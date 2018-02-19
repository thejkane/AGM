#ifndef __EAGM_RUNTIME__
#define __EAGM_RUNTIME__

#include <boost/graph/agm/util/eagm_buckets.hpp>

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
         typename OwnerMap,
         typename message_generator,
         typename eagm_config_t,
         typename rt_wrapper_gen>
class eagm_runtime_wrapper;
      
template<typename MessageGenerator,
         typename OwnerMap>
struct eagm_runtime_wrapper_gen {
public:
  amplusplus::transport& transport;
  MessageGenerator& msg_gen;
  instance_params& rtparams;
  OwnerMap& owner;
  int core_offset;
  
  eagm_runtime_wrapper_gen(amplusplus::transport& _t,
                           MessageGenerator& _mgen,
                           instance_params& _rtparams,
                           OwnerMap& _owner,
                           int offs): transport(_t),
                                      msg_gen(_mgen),
                                      rtparams(_rtparams),
                                      owner(_owner),
                                      core_offset(offs){}
                                             
  
  template <typename WorkItem,
            typename EAGMConfig,
            typename ProcessingFunction>
  struct inner {
    typedef eagm_runtime_wrapper<WorkItem,
                                 ProcessingFunction,
                                 OwnerMap,
                                 MessageGenerator,
                                 EAGMConfig,
                                 eagm_runtime_wrapper_gen<MessageGenerator,
                                                          OwnerMap> > type;
  };
};

      
template<typename work_item,
         typename processing_fn,
         typename OwnerMap,
         typename message_generator,
         typename eagm_config_t,
         typename rt_wrapper_gen>
class eagm_runtime_wrapper {
  
  // processing function
  struct handler_function;
  typedef eagm_runtime_wrapper<work_item,
                               processing_fn,
                               OwnerMap,
                               message_generator,
                               rt_wrapper_gen,
                               eagm_config_t> self_type;
  
public:
  typedef typename message_generator::template call_result<work_item, 
                                                           handler_function, 
                                                           work_item_owner<OwnerMap, work_item>,
                                                           amplusplus::no_reduction_t>::type RelaxMessage;

  // all the buckets
  typedef eagm_buckets<work_item,
                       typename eagm_config_t::global_ordering_t,
                       eagm_config_t,
                       RelaxMessage> buckets_t;


  
  eagm_runtime_wrapper(processing_fn& _pf,
                       eagm_config_t& _eagmconfig,
                       rt_wrapper_gen& _gen):
    dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item>(), 0)),
    transport(_gen.transport),
    rt_params(_gen.rtparams),
    pf(_pf),
    eagm_config(_eagmconfig),
    core_offset(_gen.core_offset),
    owner(_gen.owner),
    should_break(false),
    relax_msg(_gen.msg_gen,
              transport,
              work_item_owner<OwnerMap, work_item>(owner),
              amplusplus::no_reduction),
    all_buckets(eagm_config.global_ord,
                eagm_config,
                boost::graph::agm::spatial_level::global,
                transport,
                relax_msg)    
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

    //    debug("Decreasing activity count by 1");
    // done processing the workitem
    transport.decrease_activity_count(1);
    //debug("[[[[[Decreased activity count by 1]]]]]");
  }

  template<typename WIIterator>
  void push_work_items(WIIterator _start,
                       WIIterator _end) {
    all_buckets.push(_start, _end);
  }

  void push_work_item(work_item& wi,
                      int thread) {
    all_buckets.push(wi, thread);
  }

  bool all_buckets_empty() {
    return all_buckets.empty();
  }


  void process_bucket(int tid) {
    info("inside process bucket ! ");
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

      debug("Calling begin epoch");
      transport.begin_epoch();

      debug("Done calling begin epoch");
      for (typename std::vector<work_item>::size_type i = tid ;
           i < initial_workitems.size(); i+= nthreads) {
        const work_item& wi = initial_workitems[i];
        transport.increase_activity_count(1);        
        pf(wi, tid, all_buckets);
        transport.decrease_activity_count(1);
      }

      debug("Done pushing initial work items, entering to danger zone...");
      t_bar->wait();
      
      // execution starts
      time_type start = get_time();

      all_buckets.process(tid);

      time_type end = get_time();
      elapsed_time = (end-start);

      debug("$$$$$$$$$$$$$$$$ ENDING $$$$$$$$$$$$$$$$$$$$$$");
    }
  }


private:
  const int dummy_first_member_for_init_order; // Unused
  amplusplus::transport transport;
  instance_params rt_params;
  processing_fn& pf;
  eagm_config_t& eagm_config;
  int core_offset;
  OwnerMap& owner;
  volatile bool should_break;
  RelaxMessage relax_msg;
  buckets_t all_buckets;  
  shared_ptr<amplusplus::detail::barrier> t_bar;
  time_type elapsed_time;
};

#define EAGM_RUNTIME_WRAPPER_PARAMS \
      typename work_item, typename ordering, typename processing_fn, typename OwnerMap, typename message_generator, typename rt_wrapper_gen

#define EAGM_RUNTIME_WRAPPER_TYPE \
      eagm_runtime_wrapper<work_item, ordering, processing_fn, OwnerMap, message_generator, rt_wrapper_gen>


template<EAGM_RUNTIME_WRAPPER_PARAMS>
struct EAGM_RUNTIME_WRAPPER_TYPE::
handler_function {

  handler_function() : self(NULL){}
  handler_function(eagm_runtime_wrapper& _self) : self(&_self) {}

  void operator() (const work_item& wi) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->receive(wi, tid);
  }
  
protected:
  eagm_runtime_wrapper* self;
};

      
}}}
#endif
