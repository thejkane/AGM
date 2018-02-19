#ifndef __AGM_RUNTIME__
#define __AGM_RUNTIME__

#include <boost/graph/agm/util/bucket.hpp>

#include "runtime_common.hpp"

namespace boost { namespace graph { namespace agm {

      // TODO We need to incorporate a partition finder
      // and a bucket structure needs to be templated.
      // Finding buckets by the strict weak ordering
      // relation is inefficient due to synchronization
      // and locks.
      // PartitionFinder -- Given a work-item generated
      // using a work item in 0th bucket, this function
      // will say the the bucket id that work item should gox

// generator
template<typename work_item,
         typename ordering,
         typename processing_fn,
         typename OwnerMap,
         typename message_generator,
         typename rt_wrapper_gen>
class runtime_wrapper;
      
template<typename MessageGenerator,
         typename OwnerMap>
struct runtime_wrapper_gen {
public:
  amplusplus::transport& transport;
  MessageGenerator& msg_gen;
  instance_params& rtparams;
  OwnerMap& owner;
  int core_offset;
  
  runtime_wrapper_gen(amplusplus::transport& _t,
                      MessageGenerator& _mgen,
                      instance_params& _rtparams,
                      OwnerMap& _owner,
                      int offs): transport(_t),
                                 msg_gen(_mgen),
                                 rtparams(_rtparams),
                                 owner(_owner),
                                 core_offset(offs) {}
  
  template <typename WorkItem,
            typename Ordering,
            typename ProcessingFunction>
  struct inner {
    typedef runtime_wrapper<WorkItem,
                            Ordering,
                            ProcessingFunction,
                            OwnerMap,
                            MessageGenerator,
                            runtime_wrapper_gen<MessageGenerator, OwnerMap> > type;
  };
};

      
template<typename work_item,
         typename ordering,
         typename processing_fn,
         typename OwnerMap,
         typename message_generator,
         typename rt_wrapper_gen>
class runtime_wrapper {
  
  // processing function
  struct handler_function;
  typedef runtime_wrapper<work_item,
                          ordering,
                          processing_fn,
                          OwnerMap,
                          message_generator,
                          rt_wrapper_gen> self_type;
  
public:
  // all the buckets
  typedef buckets<work_item, ordering> buckets_t;
  // single bucket
  typedef typename bucket_structure<work_item>::Bucket a_bucket_t;
  typedef typename message_generator::template call_result<work_item, 
                                                           handler_function, 
                                                           work_item_owner<OwnerMap, work_item>,
                                                           amplusplus::no_reduction_t>::type RelaxMessage;

  
  runtime_wrapper(processing_fn& _pf,
                  ordering& _ord,
                  rt_wrapper_gen& _gen):
    dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item>(), 0)),
    transport(_gen.transport),
    rt_params(_gen.rtparams),
    pf(_pf),
    core_offset(_gen.core_offset),
    owner(_gen.owner),
    should_break(false),
    all_buckets(transport.size(), _ord),
    relax_msg(_gen.msg_gen,
              transport,
              work_item_owner<OwnerMap, work_item>(owner),
              amplusplus::no_reduction)
  {
    initialize();
  }
  
  void initialize() {
    relax_msg.set_handler(handler_function(*this));
  }

  void receive(const work_item& wi, int tid) {
    //call the processing function
    // put out set to the bucket
    //std::vector<work_item> out;
    debug("receving.");
    pf(wi, tid, all_buckets);
  }

  template<typename WIIterator>
  void push_work_items(WIIterator _start,
                       WIIterator _end) {
    all_buckets.push(_start, _end);
  }

  void push_work_item(work_item& wi) {
    all_buckets.push(wi);
  }

  bool all_buckets_empty() {
    return all_buckets.empty();
  }

  a_bucket_t*
  get_next_bucket() {
    if (all_buckets.empty()) {
      return NULL;
    } else
      return all_buckets.get_current_bucket();
  }

  void process_bucket(a_bucket_t* _bkt,
                      int tid) {

    int nthreads = rt_params.threads;
    // bucket has parallel work
    typename a_bucket_t::size_type current_bucket_start, current_bucket_end;
    current_bucket_start = 0;
    current_bucket_end = _bkt->size();

    do {
      unsigned long all_starting_sizes;
      const unsigned long starting_size = current_bucket_end - current_bucket_start;
      t_bar->wait();
      
      {
        amplusplus::scoped_epoch_value epoch(transport, starting_size, all_starting_sizes);
      
        while(current_bucket_start != current_bucket_end) {
          for (typename a_bucket_t::size_type i = current_bucket_start + tid ;
               i < current_bucket_end ; i+= nthreads) {
            work_item& wi = (*_bkt)[i];
            relax_msg.send(wi);
          }

	  t_bar->wait();
          current_bucket_start = current_bucket_end;
          current_bucket_end = _bkt->size();
	  t_bar->wait();
          // need a thread barrier ?
          
          info("Bucket start : " , current_bucket_start);
          info("Bucket end : " , current_bucket_end);
        }
      }

      if (all_starting_sizes == 0) {
        // all the processes are done with the current bucket
        if (tid == 0) {
          // done processing the  bucket
          _bkt->clear();
          all_buckets.pop_current_bucket();
        }
        
        t_bar->wait();

        break;
      } else
        current_bucket_end = _bkt->size();

    } while(true);
  }

  // returns true if there are no more buckets
  // to process
  bool synchronize(){
    return all_buckets.synchronize_current_bucket();
  }

  time_type get_elapsed_time() {
    return elapsed_time;
  }

  void set_threads(int t) {
    transport.set_nthreads(t);
  }

  void operator()(int tid,
                  std::vector<work_item>& initial_workitems) {
    int nthreads = rt_params.threads;

    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      if (tid == 0)
        t_bar.reset(new amplusplus::detail::barrier(nthreads));

      // This barrier acts as a temporary barrier until we can be sure t_bar is initialized 
      { amplusplus::scoped_epoch epoch(transport); } 
      // Now above if branch needs to be executed to every thread
      // Therefore wait till every thread comes to this point
      t_bar->wait();

      // if two processes are running on the same node, core_offset
      // is important to achieve thread affinity
      if (pin(tid+core_offset) != 0) {
        std::cerr << "[ERROR] Unable to pin current thread to "
                  << "core : " << tid << std::endl;
        assert(false);
      }

      // wait till all threads are pinned
      t_bar->wait();
      { amplusplus::scoped_epoch epoch(transport); }

      validate_thread_core_relation();

      for (typename std::vector<work_item>::size_type i = tid ;
           i < initial_workitems.size(); i+= nthreads) {
        work_item& wi = initial_workitems[i];
        push_work_item(wi);
      }

      debug("Done pushing initial work items, entering to danger zone...");
      t_bar->wait();
      
      // execution starts
      time_type start = get_time();
      
      while(true) {
        a_bucket_t* top_bucket = get_next_bucket();
        if (top_bucket != NULL) {
          process_bucket(top_bucket, tid);
        } else {
          debug("Top bucket is null ", tid);
        }

        if (tid == 0) {
          if (synchronize())
            should_break = true;
        }
        t_bar->wait();

        if (should_break)
          break;
      }

      time_type end = get_time();
      elapsed_time = (end-start);
    }
  }


private:
  const int dummy_first_member_for_init_order; // Unused
  amplusplus::transport transport;
  instance_params rt_params;
  processing_fn& pf;
  int core_offset;
  OwnerMap& owner;
  volatile bool should_break;
  buckets_t all_buckets;
  RelaxMessage relax_msg;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  time_type elapsed_time;
};

#define RUNTIME_WRAPPER_PARAMS \
      typename work_item, typename ordering, typename processing_fn, typename OwnerMap, typename message_generator, typename rt_wrapper_gen

#define RUNTIME_WRAPPER_TYPE \
      runtime_wrapper<work_item, ordering, processing_fn, OwnerMap, message_generator, rt_wrapper_gen>


template<RUNTIME_WRAPPER_PARAMS>
struct RUNTIME_WRAPPER_TYPE::
handler_function {

  handler_function() : self(NULL){}
  handler_function(runtime_wrapper& _self) : self(&_self) {}

  void operator() (const work_item& wi) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->receive(wi, tid);
  }
  
protected:
  runtime_wrapper* self;
};

      
}}}
#endif
