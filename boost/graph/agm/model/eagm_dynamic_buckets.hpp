#ifndef __EAGM_BUCKET_UTILS__
#define __EAGM_BUCKET_UTILS__

#include <string>
#include <iostream>
#include <list>
#include <unordered_set>
#include <random>
#include <tuple>
#include <type_traits>
#include <assert.h>
#include <algorithm>
#include <atomic>

#include <boost/thread/locks.hpp>
#include <boost/parallel/append_buffer.hpp>

#include <am++/make_mpi_datatype.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/util/eagm_config.hpp>

#include "bucket_stats.hpp"

      
template<typename work_item,
         typename strict_weak_ordering,
         typename eagm_configs_t,
         typename RelaxMessage>
class eagm_buckets;      
      
template<typename work_item,
         typename eagm_configs_t,
         typename RelaxMessage>
class eagm_bucket_structure {
  /**
   * The representation : This work item
   * is used to do the comparison with the strict weak
   * ordering relation. A bucket structure must be created by providing
   * a representation.
   **/

private:

  void process_global(int tid);

  void process_node(int tid);

  void process_pqs(int tid,
                   spatial_level processing_level);

  
public:

  typedef append_buffer<work_item, 10u> Bucket;

  typedef std::priority_queue<work_item, std::vector<work_item>,
                              typename eagm_configs_t::thread_ordering_t> DefaultPriorityQueueType;
  
  typedef eagm_buckets<work_item,
                       typename eagm_configs_t::node_ordering_t,
                       eagm_configs_t,
                       RelaxMessage> node_buckets_t;

  typedef std::vector< eagm_buckets<work_item,
                                    typename eagm_configs_t::numa_ordering_t,
                                    eagm_configs_t,
                                    RelaxMessage>* > numa_buckets_t;

  typedef std::vector< DefaultPriorityQueueType* > thread_buckets_t;

  eagm_bucket_structure(work_item& rep,
                        spatial_level _level,
                        amplusplus::detail::barrier* _tbar,
                        std::vector<amplusplus::detail::barrier*>& _numabar,
                        eagm_configs_t& _config,
                        amplusplus::transport& _trans,
                        RelaxMessage& _relax) : representation(rep),
                                                current_level(_level),
                                                t_bar(_tbar),
                                                numa_barriers(_numabar),
                                                eagm_config(_config),
                                                transport(_trans),
                                                relax_msg(_relax),
                                                p_node_buckets(NULL){
    initialize();
  }

  void initialize();  
  bool empty();
  void clear();
  void push_back(work_item& wi, int tid);
  void print();

  void process(int tid,
               spatial_level processing_level);
  
  int find_numa_id(int tid) {
    // TODO
    assert(false);
  }

  ~eagm_bucket_structure() {
    clear();
  }
  
  const work_item& get_representation() {
    return representation;
  }

  Bucket* get_bucket() {
    return buffer.get();
  }

private:
  boost::shared_ptr<Bucket> buffer;
  work_item representation;
  spatial_level current_level;
  node_buckets_t* p_node_buckets;
  numa_buckets_t numa_buckets;
  thread_buckets_t thread_buckets;
  amplusplus::detail::barrier* t_bar;
  std::vector<amplusplus::detail::barrier*>& numa_barriers;
  eagm_configs_t& eagm_config;
  amplusplus::transport& transport;
  RelaxMessage& relax_msg;  

};

#define EAGM_BUCKET_STRUCTURE_PARMS \
      typename work_item, typename eagm_configs_t, typename RelaxMessage

#define EAGM_BUCKET_STRUCTURE_TYPE \
      eagm_bucket_structure<work_item, eagm_configs_t, RelaxMessage>
      
      
template<typename work_item,
         typename strict_weak_ordering,
         typename eagm_configs_t,
         typename RelaxMessage>
class eagm_buckets {

  // concept assert on whether work_item has < operators defined
public:
  typedef eagm_bucket_structure < work_item, eagm_configs_t, RelaxMessage >  bucket_ds_t;

private:
  std::list< bucket_ds_t* > all_buckets;
  uint64_t eagm_index; // only incrementing
  strict_weak_ordering swo;
  eagm_configs_t& eagm_config;
  spatial_level level;
  MPI_Datatype dt;
  typename std::list< bucket_ds_t* >::iterator current_bucket;
  bucket_stats stats;
  boost::shared_mutex rwmtx;
  amplusplus::detail::barrier* t_bar;
  std::vector< amplusplus::detail::barrier* > numa_barriers;
  volatile bool should_break;
  amplusplus::transport& transport;
  RelaxMessage& relax_msg;  
  
  void initialize() {
    current_bucket = all_buckets.begin();

    //TODO thread barriers and numabarriers should be passed from an upper class
    t_bar = new amplusplus::detail::barrier(eagm_config.threads);
    // TODO : initialize numa barriers
    numa_barriers.resize(eagm_config.numa_domains);
  }
  
public:
  eagm_buckets(strict_weak_ordering& _swo,
               eagm_configs_t& _config,          
               spatial_level _level,
               amplusplus::transport& _trans,
               RelaxMessage& _relax) : eagm_index(0),
                                       swo(_swo),
                                       eagm_config(_config),
                                       level(_level),
                                       should_break(false),
                                       transport(_trans),
                                       relax_msg(_relax),
                                       dt(amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item>())){// TODO : Must pass MPI data type from top !! 
    debug("spatial level in constructor:", level);
    initialize();
  }

  
  ~eagm_buckets() {
    clear();

    // delete barriers
    delete t_bar;
    t_bar = NULL;

    for (int i=0; i < numa_barriers.size(); ++i) {
      delete numa_barriers[i];
    }
  }

  void clear() {
    typename std::list< bucket_ds_t* >::iterator ite = all_buckets.begin();
    while (ite != all_buckets.end()) {
      (*ite)->clear();
      delete (*ite);
      ++ite;
    }

    all_buckets.clear();

  }

  void print_stats() {
    stats.print();
  }

  bool empty() {
    return (all_buckets.size() == 0);
  }

  typename bucket_ds_t::Bucket*
  get_current_bucket() {
    assert(current_bucket != all_buckets.end());
    debug("Inside get_current_bucket ...");
    return (*current_bucket)->get_bucket();
  }

  // No other threads should work while
  // popping the current bucket
  // therefore, we dont need to lock this
  void pop_current_bucket(bool gotobegin = false) {
    assert(all_buckets.size() != 0);
    (*current_bucket)->clear();
    delete (*current_bucket);
    current_bucket = all_buckets.erase(current_bucket);
    // if current_bucket points to end, then point it
    // to the begining
    if (gotobegin ||
        (current_bucket == all_buckets.end()))
      current_bucket = all_buckets.begin();
  }

  void process(int tid) {

    // if this is global level, get the current
    // bucket and process it
    if (level == spatial_level::global) {
      spatial_level child_level = eagm_config.get_child_spatial_level(level);
      if (child_level == spatial_level::nil) { // this way is more efficient
        while(true) {        
          (*current_bucket)->process(tid, level);

          if (tid == 0) {
            // we come here because current bucket is empty
            pop_current_bucket();

            if (synchronize_current_bucket()) {
              should_break = true;
            }
          }
          t_bar->wait();

          if (should_break)
            break;
        }
      } else {

        debug("Global but child level not nil");
        amplusplus::transport::end_epoch_request request = transport.i_end_epoch();
        // let everything flow asynchronounsly
        while(true) {

          while(true) {
            debug("second while loop");
            (*current_bucket)->process(tid, level);

            //debug("calling i_end_epoch");
            //transport.i_end_epoch();
            debug("called i_end_epoch");
            
            if (request.test()) {
              debug("Breaking from the loop....");
              break;
            }
          }

          if (tid == 0) {
            // we come here because current bucket is empty
            pop_current_bucket();

            if (synchronize_current_bucket()) {
              should_break = true;
            }
          }
          t_bar->wait();

          if (should_break)
            break;
        }
        
      }
    }
  }

  // returns true if current bucket is empty in
  // all ranks
  bool synchronize_current_bucket() {

    //if the global ordering is chaotic we do not
    // need to do the synchronization
    //    static_if(std::is_same<CHAOTIC_ORDERING_T, strict_weak_ordering>::value){
    
    // is current bucket in local node null ? if not get
    // the representation
    // do a MPI all to all
    std::vector<work_item> recvwis;
    recvwis.resize(eagm_config.ranks);
    work_item tosend;
    
    if ((current_bucket != all_buckets.end()) &&
        (!(*current_bucket)->empty())) {
      tosend = (*current_bucket)->get_representation();
      // MPI all gather api is confusing !
      // The receive buffer size should indicate the number of
      // elements to be received per rank not the total buffer
      // size !
    } else {
      // Assuming work item first element is a vertex
      std::get<0>(tosend) = UINT64_MAX;
    }

    MPI_Allgather((void*)&tosend, 1, dt, &recvwis[0], 1, dt, MPI_COMM_WORLD);        

    bool all_empty = true;
    for (int i=0; i < eagm_config.ranks; ++i) {
      if (std::get<0>(recvwis[i]) != UINT64_MAX) {
        empty_push(recvwis[i]);
        all_empty = false;
      }     
    }
    return all_empty;
  }


  // For debugging
#ifdef PRINT_DEBUG
  void print_wi(pwi w) {
    std::cout << " push v : " << std::get<0>(w) << std::endl;
  }

  void print_wi(int w) {
    std::cout << " push v : " << w << std::endl;
  }
#endif

  void push(work_item& wi, int tid) {

    //    debug("Increasing activity count by 1");
    transport.increase_activity_count(1);
    
    while(true) {
      typename std::list< bucket_ds_t* >::iterator ite, it;
      typename std::list< bucket_ds_t* >::size_type count, step;
      uint64_t localindex = eagm_index;
      
      // Reads
      {
        boost::shared_lock<boost::shared_mutex> lock(rwmtx);

        // find the lower bound (this is lograthmic compared to linear search)
        count = all_buckets.size();
        ite = all_buckets.begin();        
        while(count > 0) {
          it = ite;
          step = (count / 2);
          std::advance(it, step);
          const work_item& representation = (*it)->get_representation();
          if (swo(representation, wi)) {
            ite = ++it;
            count -= step + 1;
          } else
            count = step;          
        }

        // if ite is not the end then it is a lowe bound
        if (ite != all_buckets.end()) {
          const work_item& representation = (*ite)->get_representation();
          if (!swo(wi, representation) && !swo(representation, wi)) {
            // not comparable insert to the current bucket
            (*ite)->push_back(wi, tid);
            return;
          } else {
            // this means there is no bucket for wi
            localindex = eagm_index;
          }
        }
      }

      // Writes
      {
        boost::unique_lock< boost::shared_mutex > lock(rwmtx);        

        if (localindex != eagm_index)
          continue;
        
        bool wasempty = false;
        if (all_buckets.size() == 0)
          wasempty = true;
        
        // value change successful
        bucket_ds_t* bucket = new bucket_ds_t(wi,
                                              eagm_config.get_child_spatial_level(level),
                                              t_bar,
                                              numa_barriers,
                                              eagm_config,
                                              transport,
                                              relax_msg);
        bucket->push_back(wi, tid);
        all_buckets.insert(ite, bucket);

        // should be inside an ifdef
        stats.increment_buckets();
        // should be inside an ifdef
        
        if (wasempty)
          current_bucket = all_buckets.begin();

        ++eagm_index;
        return;
      }
    }
  }

  //Note : This should only be used while synchronization
  // If we are synchronizing buckets globally
  // then every current_bucket should point to the
  // begining of the bucket list.
  // TODO this needs to be changed !
  void empty_push(work_item& wi) {

    typename std::list< bucket_ds_t* >::iterator ite;
    // first see whether there is a bucket for wi,
    // if so push that to the bucket.
    ite = current_bucket;
    while (ite != all_buckets.end()) {
      bucket_ds_t* pbucket = (*ite);
      // bucket can be empty
      const work_item& representation = pbucket->get_representation();
      // is wi and representation are comparable ?
      if (!swo(wi, representation) && !swo(representation, wi)) {
        return;
      } else if (swo(representation, wi)) {
        // if representation < wi, then move to next bucket
        ++ite;
      } else {
        // this means wi < representation, then break
        break;
      }
    }

    bool wasempty = false;
    if (all_buckets.size() == 0)
      wasempty = true;
        
    // value change successful
    bucket_ds_t* bucket = new bucket_ds_t(wi,
                                          eagm_config.get_child_spatial_level(level),
                                          t_bar,
                                          numa_barriers,
                                          eagm_config,
                                          transport,
                                          relax_msg);
    all_buckets.insert(ite, bucket);

    // should be inside an ifdef
    stats.increment_buckets();
    // should be inside an ifdef
        
    if (wasempty)
      current_bucket = all_buckets.begin();
  }

  
  template <typename WorkItemIterator>
  void
  push(WorkItemIterator begin, WorkItemIterator end) {
    while(begin != end) {
      push(*begin);
      ++begin;
    }
  }

  template<typename a_bucket_t>
  void print_bucket_elements(a_bucket_t* _bkt) {
    typename a_bucket_t::size_type current_bucket_start, current_bucket_end;
    current_bucket_start = 0;
    current_bucket_end = _bkt->size();

    for (typename a_bucket_t::size_type i = current_bucket_start;
         i < current_bucket_end ; i++) {
      work_item& wi = (*_bkt)[i];
      std::cout << "(" << std::get<0>(wi) << ", " << std::get<1>(wi) << "), ";
    }
  }

  uint64_t get_buckets_created() {
    return stats.get_number_of_bkst_created();
  }

  void print() {
    stats.print();
    typename std::list< bucket_ds_t* >::iterator ite = all_buckets.begin();
    int j=0;
    for (; ite != all_buckets.end(); ++ite) {
      bucket_ds_t* pbucket = (*ite);
      typename bucket_ds_t::Bucket* theBucket = pbucket->get_bucket();
      std::cout << "Bucket : " << j
                << " [REP : " << "(" << std::get<0>(pbucket->get_representation())
                << ", " << std::get<1>(pbucket->get_representation()) << ")]"
                << ", elements : " << std::endl;
      //print_bucket_elements(theBucket);
      pbucket->print();
      std::cout << std::endl;

      ++j;
    }
  }
};



template<EAGM_BUCKET_STRUCTURE_PARMS>  
void
EAGM_BUCKET_STRUCTURE_TYPE::initialize() {
  
  if (current_level == node) {
    // get ordering type for node
    p_node_buckets = new eagm_buckets<work_item,
                                      typename eagm_configs_t::node_ordering_t,
                                      eagm_configs_t,
                                      RelaxMessage>(eagm_config.node_ord,
                                                    eagm_config,
                                                    current_level,
                                                    transport,
                                                    relax_msg);
  } else if (current_level == numa) {
      
    numa_buckets.resize(eagm_config.numa_domains);

    for(int i=0; i < eagm_config.numa_domains; ++i) {
      // get ordering type for numa
      numa_buckets[i] = new eagm_buckets<work_item,
                                         typename eagm_configs_t::numa_ordering_t,
                                         eagm_configs_t,
                                         RelaxMessage>(eagm_config.numa_ord,
                                                       eagm_config,
                                                       current_level,
                                                       transport,
                                                       relax_msg);
    }

  } else if (current_level == thread) {

    thread_buckets.resize(eagm_config.threads);

    for(int i=0; i < eagm_config.threads; ++i) {
      thread_buckets[i] = new DefaultPriorityQueueType(eagm_config.thread_ord);
    }
    
  } else if (current_level == nil) {

    boost::shared_ptr<Bucket> p(new Bucket);
    buffer.swap(p);
    
  } else {
    error("Unknown spatial level", current_level);
    assert(false);
  }
}

  
template<EAGM_BUCKET_STRUCTURE_PARMS>
void
EAGM_BUCKET_STRUCTURE_TYPE::process_global(int tid) {
  int nthreads = eagm_config.threads;
  // bucket has parallel work
  typename Bucket::size_type current_bucket_start, current_bucket_end;
  current_bucket_start = 0;
  current_bucket_end = buffer->size();

  do {
    unsigned long all_starting_sizes;
    const unsigned long starting_size = current_bucket_end - current_bucket_start;
    t_bar->wait();
      
    {
      amplusplus::scoped_epoch_value epoch(transport, starting_size, all_starting_sizes);
      
      while(current_bucket_start != current_bucket_end) {
        for (typename Bucket::size_type i = current_bucket_start + tid ;
             i < current_bucket_end ; i+= nthreads) {
          work_item& wi = (*buffer)[i];
          relax_msg.send(wi);
        }

        t_bar->wait();
        current_bucket_start = current_bucket_end;
        current_bucket_end = buffer->size();
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
        buffer->clear();
      }
        
      t_bar->wait();

      break;
    } else
      current_bucket_end = buffer->size();

  } while(true);
}


template<EAGM_BUCKET_STRUCTURE_PARMS>
void
EAGM_BUCKET_STRUCTURE_TYPE::process_node(int tid) {
  // process the buffer
  int nthreads = eagm_config.threads;
  // bucket has parallel work
  typename Bucket::size_type current_bucket_start, current_bucket_end;
  current_bucket_start = 0;
  current_bucket_end = buffer->size();
  t_bar->wait();
    
  while(current_bucket_start != current_bucket_end) {
    for (typename Bucket::size_type i = current_bucket_start + tid ;
         i < current_bucket_end ; i+= nthreads) {
      work_item& wi = (*buffer)[i];
      relax_msg.send(wi);
    }

    t_bar->wait();
    current_bucket_start = current_bucket_end;
    current_bucket_end = buffer->size();
    t_bar->wait();
    // need a thread barrier ?
          
    info("Bucket start : " , current_bucket_start);
    info("Bucket end : " , current_bucket_end);
  }
}

template<EAGM_BUCKET_STRUCTURE_PARMS>
void
EAGM_BUCKET_STRUCTURE_TYPE::process_pqs(int tid,
                                        spatial_level processing_level) {
  while(!thread_buckets[tid]->empty()) {
    const work_item& wi = thread_buckets[tid]->top();
    thread_buckets[tid]->pop();
    
    bool wasempty = thread_buckets[tid]->empty();

    relax_msg.send(wi);
  }
}  
  // ThreadDomainBarrier will be a barrier only for the given
  // spatial level. For example for a give numa domain
  // it will be a barrier for the threads in that numa domain
  // e.g., if threads 3 and 4  are in numa node 2
  // and if the buffer we are executing belongs to a numa
  // domain then the ThreadDomainBarrier passed will
  // only be a barrier for threads 3 and 4
  // Lets think about NUMA later !!
template<EAGM_BUCKET_STRUCTURE_PARMS>
void
EAGM_BUCKET_STRUCTURE_TYPE::process(int tid,
                                    spatial_level processing_level) {
  if (current_level == nil) {
    // there are three possibilities here -- global, node or numa
    if (processing_level == global) {
      process_global(tid);
    } else if (processing_level == node) {
      process_node(tid);
    } else if (processing_level == numa) {
      // TODO
    } else if (processing_level == thread) {
      assert(false);
    }
  } else if (current_level == node) {
    // processing level must be global
    assert(processing_level == global);
    assert(p_node_buckets != NULL);
    p_node_buckets->process(tid);
  } else if (current_level == numa) {
    // processing level must be either global or node
    assert((processing_level == global) ||
           (processing_level == node));
    for (int i=0; i < numa_buckets.size(); ++i) {
      numa_buckets[i]->process(tid);
    }
  } else if (current_level == thread) {
    debug("Current level thread");
    // thread barrier changes depending upon the processing
    // level. if global or node t_bar and if it is numa,
    // shuld use the barrier relevant to current tid
    process_pqs(tid, processing_level); // TODO
  } else {
    error("invalid current level : " , current_level);
    assert(false);
  }
}  
  
template<EAGM_BUCKET_STRUCTURE_PARMS>  
bool
EAGM_BUCKET_STRUCTURE_TYPE::empty() {

  if ((current_level == node) &&
      (p_node_buckets != NULL)) {
    return p_node_buckets->empty();
  } else if ((current_level == numa) &&
             !numa_buckets.empty()) {
    bool empty = true;
    for (int i=0; i < numa_buckets.size(); ++i) {
      empty = empty & numa_buckets[i]->empty();
      if (!empty)
        return false;
    }

    return empty;
  } else if ((current_level == thread) &&
             !thread_buckets.empty()) {
    bool empty = true;
    for (int i=0; i < thread_buckets.size(); ++i) {
      empty = empty & thread_buckets[i]->empty();
      if (!empty)
        return false;
    }

    return empty;
  } else {
    assert(current_level==nil);
    return (buffer->size()==0);
  }
}

template<EAGM_BUCKET_STRUCTURE_PARMS>  
void
EAGM_BUCKET_STRUCTURE_TYPE::clear() {

  if ((current_level == node) &&
      (p_node_buckets != NULL)) {
    p_node_buckets->clear();
    delete p_node_buckets;
    p_node_buckets = NULL;
  } else if (current_level == numa) {
    for (int i=0; i < numa_buckets.size(); ++i) {
      numa_buckets[i]->clear();
      delete numa_buckets[i];
    }
    numa_buckets.clear();    
  } else if (current_level == thread) {
    for (int i=0; i < thread_buckets.size(); ++i) {
      delete thread_buckets[i];
    }
    thread_buckets.clear();
  } else {
    //    debug("here ..", current_level);
    //assert(current_level==nil);
    if (buffer)
      buffer->clear();
  }
}

template<EAGM_BUCKET_STRUCTURE_PARMS>  
void
EAGM_BUCKET_STRUCTURE_TYPE::push_back(work_item& wi,
                                 int tid) {
  if ((current_level == node) &&
      (p_node_buckets != NULL)) {
    p_node_buckets->push(wi, tid);
  } else if ((current_level == numa) &&
             !numa_buckets.empty()) {
    // Handle pushing to NUMA buckets
    int numaid = find_numa_id(tid);
    numa_buckets[numaid]->push(wi, tid);
  } else if ((current_level == thread) &&
             !thread_buckets.empty()) {
    thread_buckets[tid]->push(wi);
  } else {
    buffer->push_back(std::move(wi));
  }
}

template<EAGM_BUCKET_STRUCTURE_PARMS>  
void
EAGM_BUCKET_STRUCTURE_TYPE::print() {
  
  if (p_node_buckets != NULL) {
    std::cout << "[NODE_START]" << std::endl;
    p_node_buckets->print();
    std::cout << "[NODE_END]" << std::endl;
  }

  if (!numa_buckets.empty()) {
    std::cout << "[NUMA_START]" << std::endl;
    for (int i=0; i < numa_buckets.size(); ++i) {
      std::cout << "{DOMAIN:" << i << "}" << std::endl;
      numa_buckets[i]->print();
      std::cout << "{DOMAIN:" << i << "-END}" << std::endl;
    }
    std::cout << "[NUMA_END]" << std::endl;
  }

  if (!thread_buckets.empty()) {
    std::cout << "[THREAD_START] " << std::endl;
    for (int i=0; i < thread_buckets.size(); ++i) {
      std::cout << "{THREAD:" << i << "}" << std::endl;
      thread_buckets[i]->print();
      std::cout << "{THREAD:" << i << "-END}" << std::endl;
    }
    std::cout << "[THREAD_END]" << std::endl;
  }

  if (buffer) {
    std::cout << "[BS]{ ";
    typename Bucket::size_type current_bucket_start, current_bucket_end;
    current_bucket_start = 0;

    current_bucket_end = buffer->size();

    for (typename Bucket::size_type i = current_bucket_start;
         i < current_bucket_end ; i++) {
      work_item& wi = (*buffer)[i];
      std::cout << "(" << std::get<0>(wi) << ", " << std::get<1>(wi) << "), ";
    }

    std::cout << " }[BE]";
  }
}

      
}}}

#endif
