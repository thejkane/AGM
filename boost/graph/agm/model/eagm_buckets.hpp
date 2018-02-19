#ifndef __EAGM_BUCKET_STRUCT__
#define __EAGM_BUCKET_STRUCT__

#include <boost/graph/agm/model/eagm_node_defs.hpp>
#include <boost/graph/agm/util/bucket_stats.hpp>

#define NO_CHANGE_IN_PROGRESS 0            
#define CHANGE_IN_PROGRESS 1      

namespace boost { namespace graph { namespace agm {

template<typename work_item,
         typename eagm_bucket_traits_t>
class eagm_buckets {
  
public:  

  typedef typename eagm_bucket_traits_t::next_level_trait_t::spatial_t spatial_t;
  typedef typename eagm_bucket_traits_t::next_level_trait_t child_level_trait_t;
  typedef eagm_node_structure <work_item,
                               typename child_level_trait_t::order_t,
                               typename child_level_trait_t::spatial_t,
                               child_level_trait_t,
                               typename child_level_trait_t::next_level_trait_t> bucket_node_t;

  typedef typename eagm_bucket_traits_t::eagm_config eagm_config_t;
  typedef typename eagm_bucket_traits_t::runtime runtime_t;
  typedef typename eagm_bucket_traits_t::order_t strict_weak_order_t;
  typedef typename eagm_bucket_traits_t::spatial_t current_level_spatial_t;
  typedef typename std::list< bucket_node_t* > bucket_node_container_t;
  typedef typename std::list< bucket_node_t* >::iterator bucket_node_container_iterator_t;

  eagm_buckets(const eagm_config_t& _config,
               runtime_t& _rt) : config(_config),
                                 rt(_rt),
                                 eagm_index(0),
                                 should_break(false),
                                 list_status(NO_CHANGE_IN_PROGRESS),
				 forward(config.forward_scheduler)
#ifdef PBGL2_PRINT_WORK_STATS
                               ,stats(rt.get_nthreads())
#endif                                 
  {

    typename eagm_bucket_traits_t::spatial_t sp;
    set_ordering(sp);
    current_bucket = all_buckets.begin();
  }

  void process(int tid) {
    spatial_t spchild;
    //std::cout << "Calling process for " << typeid(spchild).name() << std::endl;
    //std::cout << "Current level spatial type " << typeid(current_level_spatial_t).name() << std::endl;
    current_level_spatial_t sp;
    process(tid, sp);
  }

  void process(int tid, spatial_numa_buffer_tag) {
    // in this case synchronization done only
    // for threads in the NUMA domain.
    //std::cout << "NUMA coming here ." << std::endl;
    numa_process(tid);
  }

  void process(int tid, spatial_numa_tag) {
    // in this case synchronization done only
    // for threads in the NUMA domain.
    // This is for configurations like Global --> Numa --> Thread --> Nil 
    //std::cout << "NUMA coming here (thread) ." << std::endl;
    numa_process(tid);
  }
  
  void process(int tid, spatial_numa_select_pq_or_buffer_tag) {
    // in this case synchronization done only
    // for threads in the NUMA domain.
    numa_process(tid);
  }
  
  template<typename T>
  void process(int tid, T) {
    //std::cout << "Inside normal thread execution. Should not reach here for NUMA!" << std::endl;
    // while child buckets are empty continue
    // when they are empty go to the next bucket
    if (current_bucket == all_buckets.end())
      return;

    while(current_bucket != all_buckets.end()) {
      // Leaf level buckets should make sure they process
      // all the elements currently in the bucket
      (*current_bucket)->process(tid);

      // Wait till all threads finish processing the
      // current bucket. We now know all threads
      // have done processing the current bucket, therefore,
      // we switch to the next bucket
      rt.wait_for_threads_to_reach_here(tid);
      if (rt.is_main_thread(tid)) {
        // move to next bucket
	if (forward)
	  ++current_bucket;
	else {
	  (*current_bucket)->clear();
	  all_buckets.erase(current_bucket);
	  current_bucket = all_buckets.begin();
	}
	//	fprintf(stderr, "spinning in A\n");
      }
      // To make sure that no other thread, will
      // start processing before main thread switch buckets
      rt.wait_for_threads_to_reach_here(tid);
    }

    rt.wait_for_threads_to_reach_here(tid);    
    if (rt.is_main_thread(tid))
      current_bucket = all_buckets.begin();
    rt.wait_for_threads_to_reach_here(tid);
  }

  void numa_process(int tid) {
    //std::cout << "NUMA execution ..." << std::endl;
    // while child buckets are empty continue
    // when they are empty go to the next bucket
    // synchronization localized to NUMA domain
    // that thread belongs
    if (current_bucket == all_buckets.end())
      return;

    while(current_bucket != all_buckets.end()) {
      // Leaf level buckets should make sure they process
      // all the elements currently in the bucket
      (*current_bucket)->process(tid);

      // Wait till all threads finish processing the
      // current bucket. We now know all threads
      // have done processing the current bucket, therefore,
      // we switch to the next bucket
      rt.wait_for_numa_domain_threads_to_reach_here(tid);
      if (rt.is_main_thread_in_numa_domain(tid)) {
        // move to next bucket
        ++current_bucket;
      }
      rt.wait_for_numa_domain_threads_to_reach_here(tid);
    }

    rt.wait_for_numa_domain_threads_to_reach_here(tid);
    if (rt.is_main_thread_in_numa_domain(tid))
      current_bucket = all_buckets.begin();
    rt.wait_for_numa_domain_threads_to_reach_here(tid);
  }
  
  // Note : This should only be used while synchronization
  // If we are synchronizing buckets globally
  // then every current_bucket should point to the
  // begining of the bucket list. Also, note this is called
  // in single thread.
  void balance_buckets(work_item& wi) {
    bucket_node_container_iterator_t ite, it;
    typename bucket_node_container_t::size_type count, step;
    
    // find the lower bound (this is lograthmic compared to linear search)
    count = all_buckets.size();
    ite = current_bucket;
    assert(ite == all_buckets.begin());
    
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
    
    // if ite is not the end then it is a lower bound
    if (ite != all_buckets.end()) {
      const work_item& representation = (*ite)->get_representation();
      if (!swo(wi, representation) && !swo(representation, wi)) {
        // bucket already exists, just return
        return;
      } 
    }

    bucket_node_t* bucket = new bucket_node_t(wi, config, rt);
    all_buckets.insert(ite, bucket);
    current_bucket = all_buckets.begin();
  }

  // TODO Not quite sure what empty means here.
  // Is it whether current_bucket is empty or is it that
  // aull buckets are empty
  bool empty(int tid) {
    // For the moment lets only check the current bucket
    if (current_bucket != all_buckets.end())
      return (*current_bucket)->empty(tid);
  }

  // Returns true if all the buckets are processed
  // Should be invoked in a single thread
  bool synchronize_current_bucket() {

    //    std::cout << "==== Inside Synchronization =====" << std::endl;
    // is current bucket in local node null ? if not get
    // the representation
    // do a MPI all to all
    std::vector<work_item> recvwis(rt.get_nranks());
    work_item tosend;
    
    if (current_bucket != all_buckets.end()) {
      tosend = (*current_bucket)->get_representation();
    } else {
      // Assuming work item first element is a vertex
      std::get<0>(tosend) = UINT64_MAX;
    }
    
    // MPI all gather api is confusing !
    // The receive buffer size should indicate the number of
    // elements to be received per rank not the total buffer
    // size !
    rt.do_all_gather((void*)&tosend, 1, &recvwis[0], 1);

    bool all_empty = true;
    for (int i=0; i < rt.get_nranks(); ++i) {
      if (std::get<0>(recvwis[i]) != UINT64_MAX) {
        balance_buckets(recvwis[i]);
        all_empty = false;
      }     
    }

    if (!all_empty) {
      // current bucket must be the begining bucket
      assert(current_bucket == all_buckets.begin());
    }
    
    return all_empty;
  }

  void process_global_buckets(int tid) {
    // this is where processing starts
    // should only be called on the global bucket
    // a global bucket calls process in the node relevant
    // for the current bucket

    // make sure current bucket is point to a valid place
    //    if (current_bucket != all_buckets.end()) {
    // process all global buckets
    while(true) {

      if (rt.is_main_thread(tid)) {          
	if (synchronize_current_bucket()) {
	  should_break = true;
	}
      }
      // wait till main thread reach here
      // when main thread is here we know that
      // should_break is properly assigned
      rt.wait_for_threads_to_reach_here(tid);

      if (should_break)
	break;

      if (rt.is_main_thread(tid)) {
	rt.increase_activity_count(tid,
				   (*current_bucket)->
				   reduce_thread_push_counts());
      }

      // Wait till main thread increases the
      // activity count
      rt.wait_for_threads_to_reach_here(tid);
      process_current_bucket(tid);
    }
  }

  void process_current_bucket(int tid) {    
    unsigned long all_starting_sizes = 0;
    do {
      const unsigned long starting_size = rt.get_activity_count();
      rt.wait_for_threads_to_reach_here(tid);
      
      {
        runtime_epoch_value<runtime_t> epoch(rt, starting_size, all_starting_sizes);
        (*current_bucket)->process(tid);
      }

      if (all_starting_sizes == 0) {
        // all the processes are done with the current bucket
        if (rt.is_main_thread(tid)) {
          // done processing the  bucket
          // clear data for current bucket
          (*current_bucket)->clear();
          // delete the current bucket
          delete (*current_bucket);
          // assign next bucket as the current bucket
          current_bucket = all_buckets.erase(current_bucket);

#ifdef PBGL2_PRINT_WORK_STATS
          stats.increment_buckets();
#endif          
        }
        
        rt.wait_for_threads_to_reach_here(tid);
        break;
      }
    } while(true);
  }

  void process_local_buckets() {
  }

  void set_ordering(spatial_global_tag _gt) {
    swo = config.global_ord;
  }

  void set_ordering(spatial_node_tag _gt) {
    swo = config.node_ord;
  }

  void set_ordering(spatial_numa_tag _gt) {
    swo = config.numa_ord;
  }

  void set_ordering(spatial_thread_tag _gt) {
    assert(false); // should not come here !
  }  

  void set_ordering(spatial_node_select_pq_or_buffer_tag _gt) {
    swo = config.node_ord;    
  }  

  void set_ordering(spatial_numa_select_pq_or_buffer_tag _gt) {
    swo = config.numa_ord;    
  }    

  void push(const work_item& wi, int tid) {
    current_level_spatial_t cls;
    push(wi, tid, cls);
  }

  void clear() {    
    bucket_node_container_iterator_t begin = all_buckets.begin();
    for (; begin != all_buckets.end(); ++begin) {
      (*begin)->clear();
    }

    all_buckets.clear();
  }
  
  template<typename T>
  void push(const work_item& wi, int tid, T) {
    static_assert((!std::is_same<T, spatial_global_tag>::value),
                  "This push function can only be called for non-global buckets.");
    general_push(wi, tid);
  }

  void push(const work_item& wi, int tid, spatial_global_tag) {
    // Note : This is bad, because we are duplicating code.
    // However, with this we avoid certain runtime checks
    // , cos we know this is invoked only for the global buckets
    // Leading to improve the performance.

    // We need to carefully handle termination.
    // The activity count should be increased only if we are
    // inserting an element to the currently processing global bucket.
    // Otherwise, we need to increase thread counts for the
    // future buckets. Also, this is global level, therefore,
    // there cant be buckets inserted before current_bucket iterator.
#ifdef PBGL2_PRINT_WORK_STATS
    stats.increment_pushes(tid);
#endif
    while(true) {
      
      bucket_node_container_iterator_t ite, it;
      typename bucket_node_container_t::size_type count, step;
      uint64_t localindex;
      
      // Reads
      {
        localindex = eagm_index;
        //        fprintf(stderr, "rli : %" PRIu64 ",  ei : %" PRIu64 "tid : %d\n",
        //        localindex, eagm_index, tid);
        
        //boost::shared_lock<boost::shared_mutex> lock(rwmtx);
        // std::list insert says that no iterators will
        // invalidated. Therefore, do I really need the lock ?
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

        // if ite is not the end then it is a lower bound
        if (ite != all_buckets.end()) {
          const work_item& representation = (*ite)->get_representation();
          if (!swo(wi, representation) && !swo(representation, wi)) {
            // not comparable insert to the current bucket
            //is ite same as current_bucket ? yes -- increase activity count
            // no -- increase thread_push_counts
            if (ite == current_bucket) {
              rt.increase_activity_count(tid, 1);
            } else {
              (*ite)->increment_thread_push_counts(tid);
            }
            
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
        //        fprintf(stderr, "wli : %" PRIu64 ",  ei : %" PRIu64 "tid : %d\n",
        //        localindex, eagm_index, tid);

        //boost::lock_guard<boost::mutex> lock(rwmtx);        
        int expected = NO_CHANGE_IN_PROGRESS;
        int valtoset = CHANGE_IN_PROGRESS;

        bool continue_from_begining = false;
        
        while(expected == NO_CHANGE_IN_PROGRESS) {
          if (__atomic_compare_exchange_n(&list_status,
                                          &expected,
                                          valtoset,
                                          true,
                                          __ATOMIC_RELAXED,
                                          __ATOMIC_RELAXED)) {
            // current thread got the control
            // now list_status = CHANGE_IN_PROGRESS;
            //    fprintf(stderr, "I am here ...%d\n", tid);
            assert(list_status == CHANGE_IN_PROGRESS);
            break;
          }
          
          if (localindex != eagm_index) {
            continue_from_begining = true;
            break;
          }

          expected = NO_CHANGE_IN_PROGRESS;
        }

        
        //        fprintf(stderr, "li : %" PRIu64 ",  ei : %" PRIu64 "tid : %d\n",
        //        localindex, eagm_index, tid);

        if (continue_from_begining) {          
          continue;
        }
        
        // only one thread should come here
        if (localindex != eagm_index) {
          list_status = NO_CHANGE_IN_PROGRESS;
          continue;
        }
        
        
        bool wasempty = false;
        if (all_buckets.size() == 0)
          wasempty = true;

        // value change successful
        bucket_node_t* bucket = new bucket_node_t(wi, config, rt);
        bucket->push_back(wi, tid);
        all_buckets.insert(ite, bucket);
	bucket->increment_thread_push_counts(tid);
        // TODO add stats
        
        if (wasempty) {
          current_bucket = all_buckets.begin();
	  //  rt.increase_activity_count(tid, 1);
        } /*else {
          bucket->increment_thread_push_counts(tid);
	  }*/
        
        ++eagm_index;
        list_status = NO_CHANGE_IN_PROGRESS;        
        return;
      }
    }
  }
  
  void general_push(const work_item& wi, int tid) {    
    
    while(true) {
      
      bucket_node_container_iterator_t ite, it;
      typename bucket_node_container_t::size_type count, step;
      uint64_t localindex;
      
      // Reads
      {
        localindex = eagm_index;

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

        // if ite is not the end then it is a lower bound
        if (ite != all_buckets.end()) {
          const work_item& representation = (*ite)->get_representation();
          if (!swo(wi, representation) && !swo(representation, wi)) {
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
        int expected = NO_CHANGE_IN_PROGRESS;
        int valtoset = CHANGE_IN_PROGRESS;

        bool continue_from_begining = false;
        
        while(expected == NO_CHANGE_IN_PROGRESS) {
          if (__atomic_compare_exchange_n(&list_status,
                                          &expected,
                                          valtoset,
                                          true,
                                          __ATOMIC_RELAXED,
                                          __ATOMIC_RELAXED)) {
            // current thread got the control
            // now list_status = CHANGE_IN_PROGRESS;
            //    fprintf(stderr, "I am here ...%d\n", tid);
            assert(list_status == CHANGE_IN_PROGRESS);
            break;
          }
          
          if (localindex != eagm_index) {
            continue_from_begining = true;
            break;
          }

          expected = NO_CHANGE_IN_PROGRESS;
        }

        if (continue_from_begining) {          
          continue;
        }
        
        // only one thread should come here
        if (localindex != eagm_index) {
          list_status = NO_CHANGE_IN_PROGRESS;
          continue;
        }
        
        bool wasempty = false;
        if (all_buckets.size() == 0)
          wasempty = true;
        
        // value change successful
	//	fprintf(stderr, "non-global bucket created\n");
        bucket_node_t* bucket = new bucket_node_t(wi, config, rt);
        bucket->push_back(wi, tid);
        all_buckets.insert(ite, bucket);

        // TODO add stats
        
        if (wasempty)
          current_bucket = all_buckets.begin();

        ++eagm_index;
        list_status = NO_CHANGE_IN_PROGRESS;        
        return;
      }
    }
  }

  void print() {
    bucket_node_container_iterator_t ite = all_buckets.begin();
    int j=0;
    for (; ite != all_buckets.end(); ++ite) {
      bucket_node_t* pbucket = (*ite);
      std::cout << "Bucket : " << j
                << " [REP : " << "(" << std::get<0>(pbucket->get_representation())
                << ", " << std::get<1>(pbucket->get_representation()) << ")]"
                << ", elements : " << std::endl;
      pbucket->print();
      std::cout << std::endl;

      ++j;
    }
  }

  void print_stats() {
#ifdef PBGL2_PRINT_WORK_STATS
    stats.print();
#endif    
  }
private:
  bucket_node_container_t all_buckets;
  const eagm_config_t& config;
  runtime_t& rt;
  strict_weak_order_t swo;
  uint64_t eagm_index; // only incrementing
  volatile bool should_break;
  //  boost::mutex rwmtx;
  int list_status;  
  // can list grow more than 32 bit unsigned int ?
  bucket_node_container_iterator_t current_bucket;
  bucket_node_container_iterator_t next_bucket;
  bool forward;
#ifdef PBGL2_PRINT_WORK_STATS
  bucket_stats stats;
#endif  
};

}}}

#endif
