#ifndef __EAGM_NODE_DEFS__
#define __EAGM_NODE_DEFS__

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
#include <type_traits>
#include <queue>
#include <vector>
#include <typeinfo>

#include <boost/graph/agm/model/eagm_traits.hpp>
#include <boost/graph/agm/util/priority_qs.hpp>
#include <boost/graph/agm/runtime/runtime_common.hpp>

//============================== EAGM NODE STRUCTURES ========================================
namespace boost { namespace graph { namespace agm {
      
template<typename work_item, typename current_level_trait>      
class eagm_node_base {
public:
  typedef typename current_level_trait::eagm_config eagm_config_t;
  typedef typename current_level_trait::runtime runtime_t;
  typedef typename eagm_config_t::numa_container_t numa_container_t;
  typedef typename eagm_config_t::node_container_t node_container_t;  
  
  eagm_node_base(const work_item& _rep,
                 const eagm_config_t& _config,
                 runtime_t& _rt) : representation(_rep),
                                   config(_config),
                                   rt(_rt){
    thread_push_counts.resize(rt.get_nthreads(), 0);
  }

  const work_item& get_representation() {
    return representation;
  }
  
  void increment_thread_push_counts(int tid) {
    assert(tid < rt.get_nthreads());
    ++thread_push_counts[tid];
  }

  uint64_t reduce_thread_push_counts() {
    uint64_t tot = 0;
    for(int i=0; i < rt.get_nthreads(); ++i) {
      tot += thread_push_counts[i];
    }

    return tot;
  }
  

protected:
  const work_item representation;  
  const eagm_config_t& config;
  runtime_t& rt;
  std::vector<uint64_t> thread_push_counts;  
};
      
template<typename work_item_t,
         typename ordering_t,
         typename spatial_level,         
         typename current_level_trait,
         typename next_level_trait_t>               
class eagm_node_structure {
};

// A bucket node that stores priority queues            
template<typename work_item,
         typename ordering_t,
         typename current_level_trait>               
class eagm_node_structure<work_item,
                          ordering_t,
                          boost::graph::agm::spatial_thread_tag,
                          current_level_trait,
                          boost::graph::agm::leaf_level_eagm_trait> : public eagm_node_base<work_item, current_level_trait> {

private:  
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;
  typedef std::priority_queue<work_item,
                              std::vector<work_item>, ordering_t> DefaultPriorityQueueType;
  typedef std::vector< boost::shared_ptr<DefaultPriorityQueueType> > thread_buckets_t;
  
public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                        current_level_trait>(_rep,
                                                                             _config,
                                                                             _rt){
    thread_buckets.resize(this->rt.get_nthreads());
    
    for(int i=0; i < this->rt.get_nthreads(); ++i) {
      boost::shared_ptr<DefaultPriorityQueueType>
        p(new DefaultPriorityQueueType(this->config.thread_ord));
      thread_buckets[i].swap(p);
    }
  }
  
  void push_back(const work_item& wi, int tid) {
    thread_buckets[tid]->push(wi);
  }

  void process(int tid) {
    while(!thread_buckets[tid]->empty()) {
      // TODO object copying, but tricky since we need to
      // pop also
      work_item wi = thread_buckets[tid]->top();
      thread_buckets[tid]->pop();      
      this->rt.send(wi, tid);

    }
  }

  void clear() {
    for (int i=0; i < thread_buckets.size(); ++i) {
      assert(thread_buckets[i]->empty());
    }
  }
  
  void print() {
    if (!thread_buckets.empty()) {
      std::cout << "[THREAD_START] " << std::endl;
      for (int i=0; i < thread_buckets.size(); ++i) {
        std::cout << "{THREAD:" << i << "}" << std::endl;
        print_pq(i);
        std::cout << "{THREAD:" << i << "-END}" << std::endl;
      }
      std::cout << "[THREAD_END]" << std::endl;
    }
  }

  bool empty(int tid) {
    return thread_buckets[tid]->empty();
  }

  ~eagm_node_structure() {
    clear();
  }

private:

  void print_pq(int i) {
    // disabled since this is popping elements
    /*    while(!thread_buckets[i]->empty()) {
      std::cout << "(" << std::get<0>(thread_buckets[i]->top()) << ", "
                << std::get<1>(thread_buckets[i]->top()) << "), ";
      thread_buckets[i]->pop();
    }

    std::cout << std::endl;*/
  }
  
  thread_buckets_t thread_buckets;
};

// A bucket node that sends work items in push
// represents Global (Chaotic) --> Nil
template<typename work_item,
         typename current_level_trait>               
class eagm_node_structure<work_item,
                          BufferOrdering,
                          boost::graph::agm::spatial_global_all_chaotic_buffer_tag,
                          current_level_trait,
                          boost::graph::agm::leaf_level_eagm_trait>  : public eagm_node_base<work_item, current_level_trait> {

private:  
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;

public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) : eagm_node_base<work_item,
                                                       current_level_trait>(_rep,
                                                                            _config,
                                                                            _rt){
  }
  
  void push_back(const work_item& wi, int tid) {
    this->rt.send(wi, tid);
  }

  void process(int tid) {
    // Nothing to do
  }

  void print() {
    std::cout << "[Global]{}[Global]";
  }
  
  void clear() {
  }

};


// A bucket node that stores an append buffer for a node
// represents Global --> buffer (AGM)
// This is quite close to node-->buffer implementation.
// What is the difference between node-->buffer and global-->buffer
// The process function in global-->buffer synchronizes, globally
// but synchronization in node-->buffer is local.      
template<typename work_item,
         typename current_level_trait>               
class eagm_node_structure<work_item,
                          boost::graph::agm::BufferOrdering,
                          boost::graph::agm::spatial_global_buffer_tag,
                          current_level_trait,
                          boost::graph::agm::leaf_level_eagm_trait> : public eagm_node_base<work_item, current_level_trait> {

private:  
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;

  typedef append_buffer<work_item, 10u> Bucket;
  
public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                        current_level_trait>(_rep,
                                                                             _config,
                                                                             _rt),
    current_bucket_start(0){
    
    boost::shared_ptr<Bucket> p(new Bucket);
    buffer.swap(p);
  }
  
  void push_back(const work_item& wi, int tid) {
    buffer->push_back(std::move(wi));
  }

  void process(int tid) {

    // NOTE: Clearing buffers may consume more time.
    // Therefore, think about clearing the buffer
    // after certain threashold. For example when
    // buffer reaches certain size, clear it.
    // Consider doing this we if we see bad performance
    // also measure the time for clearing.
    
    int nthreads = this->rt.get_nthreads();
    // bucket has parallel work
    typename Bucket::size_type current_bucket_end;

    // Need to make sure all threads have the same
    // current_bucket_end value. We need both barriers.
    this->rt.wait_for_threads_to_reach_here(tid);
    current_bucket_end = buffer->size();
    this->rt.wait_for_threads_to_reach_here(tid);
    
    while(current_bucket_start != current_bucket_end) {
      for (typename Bucket::size_type i = current_bucket_start + tid ;
           i < current_bucket_end ; i+= nthreads) {
        work_item& wi = (*buffer)[i];
        this->rt.send(wi, tid);
      }
      
      this->rt.wait_for_threads_to_reach_here(tid);
      if (tid == 0)
        current_bucket_start = current_bucket_end;      

      current_bucket_end = buffer->size();
      this->rt.wait_for_threads_to_reach_here(tid);

#ifdef ENABLE_INTERMEDIATE_CLEARING      
      // by here we know that all threads have the
      // same current_bucket_start and current_bucket_end
      if (current_bucket_start == current_bucket_end) {
        if (this->rt.is_main_thread(tid)) {
          buffer->clear();
        }
        
        // We need this barrier to make a thread may not go
        // and push elements to the bucket while the main
        // thread is clearing the bucket
        this->rt.wait_for_threads_to_reach_here(tid);
      }
#endif      
    }
  }

  void print() {
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

  void clear() {
    assert(buffer);
    buffer->clear();
    current_bucket_start = 0;    
  }
  
private:
  boost::shared_ptr<Bucket> buffer;
  typename Bucket::size_type current_bucket_start;
};
      

// A bucket node that stores an append buffer for a node
// represents node-->buffer      
template<typename work_item,
         typename current_level_trait>               
class eagm_node_structure<work_item,
                          boost::graph::agm::BufferOrdering,
                          boost::graph::agm::spatial_node_buffer_tag,
                          current_level_trait,
                          boost::graph::agm::leaf_level_eagm_trait> : public eagm_node_base<work_item, current_level_trait> {

private:  
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;
  typedef append_buffer<work_item, 10u> Bucket;
  
public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                        current_level_trait>(_rep,
                                                                             _config,
                                                                             _rt),
					 current_bucket_start(0){
    boost::shared_ptr<Bucket> p(new Bucket);
    buffer.swap(p);
  }
  
  void push_back(const work_item& wi, int tid) {
    buffer->push_back(std::move(wi));
  }

  void process(int tid) {
    
    int nthreads = this->rt.get_nthreads();
    // bucket has parallel work
    typename Bucket::size_type current_bucket_end;

    // Need to make sure all threads have the same
    // current_bucket_end value. We need both barriers.
    this->rt.wait_for_threads_to_reach_here(tid);
    current_bucket_end = buffer->size();
    this->rt.wait_for_threads_to_reach_here(tid);
    
    while(current_bucket_start != current_bucket_end) {
      for (typename Bucket::size_type i = current_bucket_start + tid ;
           i < current_bucket_end ; i+= nthreads) {
        work_item& wi = (*buffer)[i];
        this->rt.send(wi, tid);
      }

      //fprintf(stderr, "spinning in B\n");      
      this->rt.wait_for_threads_to_reach_here(tid);
      if (tid == 0)
	current_bucket_start = current_bucket_end;      

      current_bucket_end = buffer->size();
      this->rt.wait_for_threads_to_reach_here(tid);

#ifdef ENABLE_INTERMEDIATE_CLEARING      
      // by here we know that all threads have the
      // same current_bucket_start and current_bucket_end
      if (current_bucket_start == current_bucket_end) {
        if (this->rt.is_main_thread(tid)) {
          buffer->clear();
        }
        
        // We need this barrier to make a thread may not go
        // and push elements to the bucket while the main
        // thread is clearing the bucket
        this->rt.wait_for_threads_to_reach_here(tid);
      }
#endif 
    }
  }

  void print() {
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

  void clear() {
    buffer->clear();
    current_bucket_start = 0;    
  }

  bool empty(int tid) {
    return (buffer->size() == 0);
  }
  
private:
  boost::shared_ptr<Bucket> buffer;
  typename Bucket::size_type current_bucket_start;
};

// A bucket node that stores an append buffer for a numa domain
// represents numa-->buffer      
template<typename work_item,
         typename current_level_trait>               
class eagm_node_structure<work_item,
                          boost::graph::agm::BufferOrdering,
                          boost::graph::agm::spatial_numa_buffer_tag,
                          current_level_trait,
                          boost::graph::agm::leaf_level_eagm_trait> : public eagm_node_base<work_item, current_level_trait> {

  
private:
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;
  typedef append_buffer<work_item, 10u> Bucket;
  
public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                        current_level_trait>(_rep,
                                                                             _config,
                                                                             _rt){
    boost::shared_ptr<Bucket> p(new Bucket);
    buffer.swap(p);
  }
  
  void push_back(const work_item& wi, int tid) {
    buffer->push_back(std::move(wi));
  }

  void process(int tid) {
    int nthreads = this->rt.get_nthreads_for_numa_domain(tid);
    int index = this->rt.get_thread_index_in_numa_domain(tid);
    // bucket has parallel work
    typename Bucket::size_type current_bucket_start, current_bucket_end;
    this->rt.wait_for_numa_domain_threads_to_reach_here(tid);
    current_bucket_start = 0;
    current_bucket_end = buffer->size();
    this->rt.wait_for_numa_domain_threads_to_reach_here(tid);
    
    while(current_bucket_start != current_bucket_end) {
      for (typename Bucket::size_type i = current_bucket_start + index ;
           i < current_bucket_end ; i+= nthreads) {
        work_item& wi = (*buffer)[i];
        this->rt.send(wi, tid);
      }

      this->rt.wait_for_numa_domain_threads_to_reach_here(tid);
      current_bucket_start = current_bucket_end;
      current_bucket_end = buffer->size();
      this->rt.wait_for_numa_domain_threads_to_reach_here(tid);
      
      if (current_bucket_start == current_bucket_end) {
        if (this->rt.is_main_thread_in_numa_domain(tid)) {
          buffer->clear();
        }
      }
      this->rt.wait_for_numa_domain_threads_to_reach_here(tid);
    }
  }

  void clear() {
    buffer->clear();
  }

  bool empty(int tid) {
    return (buffer->size() == 0);
  }

  void print() {
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
  
private:
  boost::shared_ptr<Bucket> buffer;
};
      

// Decleration ...      
template<typename work_item,
         typename eagm_bucket_traits_t>
class eagm_buckets;

      
// e.g., Global --> Node --> Thread --> Nil
template<typename work_item,
         typename ordering_t,
         typename current_level_trait,
         typename next_level_trait>               
class eagm_node_structure<work_item,
                          ordering_t,
                          boost::graph::agm::spatial_node_tag,
                          current_level_trait,
                          next_level_trait> : public eagm_node_base<work_item,
                                                             current_level_trait> {

  typedef typename current_level_trait::next_level_trait_t::spatial_t nextrait_spatial;
  //nextrait must be either numa or thread
  static_assert(std::is_same<nextrait_spatial, boost::graph::agm::spatial_numa_tag>::value ||
                //             std::is_same<nextrait_spatial, boost::graph::agm::spatial_node_buffer_tag>::value ||
                std::is_same<nextrait_spatial, boost::graph::agm::spatial_numa_select_pq_or_buffer_tag>::value ||
                std::is_same<nextrait_spatial, boost::graph::agm::spatial_thread_tag>::value,
                "nextrait_spatial is not numa nor thread. something is wrong");

private:  
  // we are at the node level, create a bucket structure
  typedef eagm_buckets<work_item, current_level_trait> node_buckets_t;
  
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;


public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                        current_level_trait>(_rep,
                                                                             _config,
                                                                             _rt),
    node_buckets(_config, _rt){

  }

  
  void push_back(const work_item& wi, int tid) {
    node_buckets.push(wi, tid);
  }

  void print() {
    std::cout << "[NODE_START]" << std::endl;
    node_buckets.print();
    std::cout << "[NODE_END]" << std::endl;
  }

  void process(int tid) {
    node_buckets.process(tid);
  }

  void clear() {
    node_buckets.clear();
  }

  /*bool empty(int tid) {
    return node_buckets.empty(tid);
    }*/


private:
  node_buckets_t node_buckets;
};


// Global --> Node --> Buffer --> Nil
// Made this special, because we need to select the node level container
// type. If EAGM configuration say to use priority queues we use
// concurrent priority queues otherwise we will use append buffers.      
template<typename work_item,
         typename ordering_t,
         typename current_level_trait,
         typename next_level_trait>               
class eagm_node_structure<work_item,
                          ordering_t,
                          boost::graph::agm::spatial_node_select_pq_or_buffer_tag,
                          current_level_trait,
                          next_level_trait> : public eagm_node_base<work_item,
                                                             current_level_trait> {

  typedef typename current_level_trait::next_level_trait_t::spatial_t nextrait_spatial;
  //nextrait must be either numa or thread
  static_assert(std::is_same<nextrait_spatial, boost::graph::agm::spatial_node_buffer_tag>::value,
                "nextrait_spatial must be node buffer. something is wrong");

private:
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;
  
  // we are at the node level, create a bucket structure
  typedef eagm_buckets<work_item, current_level_trait> node_buckets_t;  
  typedef typename eagm_config_t::node_container_t node_container_t;
  typedef boost::graph::agm::concurrent_priority_queue<work_item, ordering_t> concurrent_pq_t;
  
  void initialize(buffer_container) {
    boost::shared_ptr<node_buckets_t> p(new node_buckets_t(this->config, this->rt));
    node_buckets.swap(p);
  }

  void initialize(pq_container) {
    boost::shared_ptr<concurrent_pq_t> p(new concurrent_pq_t(this->rt.get_nthreads()));
    concurrent_pq.swap(p);
  }  

  void push_back(const work_item& wi, int tid, buffer_container) {
    node_buckets->push(wi, tid);
  }

  void push_back(const work_item& wi, int tid, pq_container) {
    concurrent_pq->put(wi, tid);
  }

  void process(int tid, buffer_container) {
    node_buckets->process(tid);
  }

  void process(int tid, pq_container) {
    work_item wi;
    while(concurrent_pq->pop(wi, tid)) {
      this->rt.send(wi, tid);
    }
  }

  void clear(buffer_container) {
    node_buckets->clear();
  }

  void clear(pq_container) {
    assert(concurrent_pq->empty());
  }

  void print(buffer_container) {
    std::cout << "[NODE_START]" << std::endl;
    node_buckets->print();
    std::cout << "[NODE_END]" << std::endl;
  }

  void print(pq_container) {
    std::cout << "[NODE_START]" << std::endl;    
    /*    work_item wi;
    while(concurrent_pq->pop(wi, 0)) {
      std::cout << "(" <<  std::get<0>(wi) << ", "
                << std::get<1>(wi) << "), "; 
    }
    std::cout << std::endl;*/
    std::cout << "[NODE_END]" << std::endl;    
  }
  
public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                       current_level_trait>(_rep,
                                                                            _config,
                                                                            _rt){
    node_container_t c;
    initialize(c);
  }
  
  void push_back(const work_item& wi, int tid) {
    node_container_t c;
    push_back(wi, tid, c);
  }

  void process(int tid) {
    node_container_t c;
    process(tid, c);
  }

  void clear() {
    node_container_t c;
    clear(c);
  }

  void print() {
    node_container_t c;
    print(c);
  }


private:
  boost::shared_ptr<node_buckets_t> node_buckets;
  boost::shared_ptr<concurrent_pq_t> concurrent_pq;
};
      

// e.g., Global --> Numa --> Thread --> Nil
template<typename work_item,
         typename ordering_t,
         typename current_level_trait,
         typename next_level_trait>               
class eagm_node_structure<work_item,
                          ordering_t,
                          boost::graph::agm::spatial_numa_tag,
                          current_level_trait,
                          next_level_trait> :  public eagm_node_base<work_item,
                                                             current_level_trait> {

  typedef typename current_level_trait::next_level_trait_t::spatial_t nextrait_spatial;
  //nextrait must be either buffer or thread
  static_assert((std::is_same<nextrait_spatial, boost::graph::agm::spatial_thread_tag>::value ||
                 std::is_same<nextrait_spatial, boost::graph::agm::spatial_numa_buffer_tag>::value) ,
                "nextrait_spatial is not thread. something is wrong");

private:  
  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  // we are at the numa level, create a bucket structure for each numa domain
  typedef eagm_buckets<work_item, current_level_trait> numa_buckets_t;
  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;
  typedef std::vector<boost::shared_ptr<numa_buckets_t> > all_numa_buckets_t;

public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) :  eagm_node_base<work_item,
                                                        current_level_trait>(_rep,
                                                                             _config,
                                                                             _rt){
    numa_buckets.resize(this->rt.get_nnuma_nodes());
    for(int i=0; i < this->rt.get_nnuma_nodes(); ++i) {
      boost::shared_ptr<numa_buckets_t> p(new numa_buckets_t(this->config, this->rt));
      numa_buckets[i].swap(p); 
    }
  }
  
  void push_back(const work_item& wi, int tid) {
    // allocate memory for buckets in the same context as the tid's numa domain
    int numanode = this->rt.find_numa_node(tid);
    
    assert(numanode < numa_buckets.size());
    assert(numa_buckets[numanode]);
    
    numa_buckets[numanode]->push(wi, tid);
  }

  void print() {
    std::cout << "[NUMA_START]" << std::endl;
    for (int i=0; i < numa_buckets.size(); ++i) {
      std::cout << "{DOMAIN:" << i << "}" << std::endl;
      numa_buckets[i]->print();
      std::cout << "{DOMAIN:" << i << "-END}" << std::endl;
    }
    std::cout << "[NUMA_END]" << std::endl;
  }

  void clear() {
    // Note : does not clear the actual numa vector
    for (int i=0; i < numa_buckets.size(); ++i) {
      numa_buckets[i]->clear();
    }
  }

  void process(int tid) {
    int numanode = this->rt.find_numa_node(tid);
    numa_buckets[numanode]->process(tid);
  }

  bool empty(int tid) {
    int numanode = this->rt.find_numa_node(tid);
    numa_buckets[numanode]->empty(tid);    
  }
  
private:
  all_numa_buckets_t numa_buckets;
};


// Global --> Numa --> Buffer --> Nil
// Based on the EAGM configuration we decide whether we should
// use a concurrent priority queue for every NUMA domain or a
// separate buffer.      
template<typename work_item,
         typename ordering_t,
         typename current_level_trait,
         typename next_level_trait>               
class eagm_node_structure<work_item,
                          ordering_t,
                          boost::graph::agm::spatial_numa_select_pq_or_buffer_tag,
                          current_level_trait,
                          next_level_trait> : public eagm_node_base<work_item,
                                                             current_level_trait> {

  typedef typename current_level_trait::next_level_trait_t::spatial_t nextrait_spatial;
  //nextrait must be either buffer or thread
  static_assert(std::is_same<nextrait_spatial, boost::graph::agm::spatial_numa_buffer_tag>::value ,
                "nextrait_spatial is not numa buffer. something is wrong");

  typedef eagm_node_base<work_item,current_level_trait> base_class_t;
  
private:

  typedef typename base_class_t::eagm_config_t eagm_config_t;
  typedef typename base_class_t::runtime_t runtime_t;
  typedef typename base_class_t::numa_container_t numa_container_t;  
  
  typedef eagm_buckets<work_item, current_level_trait> numa_buckets_t;  
  typedef std::vector<boost::shared_ptr<numa_buckets_t> > all_numa_buckets_t;
  // OR
  typedef boost::graph::agm::concurrent_priority_queue<work_item, ordering_t> concurrent_pq_t;
  typedef std::vector<boost::shared_ptr<concurrent_pq_t> > all_numa_pqs_t;
  
  void initialize(buffer_container) {
    numa_buckets.resize(this->rt.get_nnuma_nodes());

    // this initialization is simple but may not be ideal. Cos
    // we do not know in which numa node the memmory will be allocated
    // for the data. But this is simple because we can avoid
    // the concurrency handling
    for (int i=0; i < this->rt.get_nnuma_nodes(); ++i) {
      boost::shared_ptr<numa_buckets_t> p(new numa_buckets_t(this->config, this->rt));
      numa_buckets[i].swap(p);
    }
  }

  void initialize(pq_container) {
    numa_pqs.resize(this->rt.get_nnuma_nodes());

    for (int i=0; i < this->rt.get_nnuma_nodes(); ++i) {
      boost::shared_ptr<concurrent_pq_t> p(new concurrent_pq_t(this->rt.get_nthreads()));
      numa_pqs[i].swap(p);
    }
  }  

  void push_back(const work_item& wi, int tid, buffer_container) {
    // allocate memory for buckets in the same context as the tid's numa domain
    int numanode = this->rt.find_numa_node(tid);
    
    assert(numanode < numa_buckets.size());
    assert(numa_buckets[numanode]);
    
    numa_buckets[numanode]->push(wi, tid);
  }

  void push_back(const work_item& wi, int tid, pq_container) {
    int numanode = this->rt.find_numa_node(tid);
    
    assert(numanode < numa_pqs.size());
    assert(numa_pqs[numanode]);

    numa_pqs[numanode]->put(wi, tid);
  }

  void process(int tid, buffer_container) {
    int numanode = this->rt.find_numa_node(tid);
    numa_buckets[numanode]->process(tid);
  }

  void process(int tid, pq_container) {
    int numanode = this->rt.find_numa_node(tid);

    while(!numa_pqs[numanode]->empty()) {
      work_item wi;
      if (numa_pqs[numanode]->pop(wi, tid)) {
        this->rt.send(wi, tid);
      }
    }
  }  
  
  void print(buffer_container) {
    std::cout << "[NUMA_START]" << std::endl;
    for (int i=0; i < numa_buckets.size(); ++i) {
      std::cout << "{DOMAIN:" << i << "}" << std::endl;
      numa_buckets[i]->print();
      std::cout << "{DOMAIN:" << i << "-END}" << std::endl;
    }
    std::cout << "[NUMA_END]" << std::endl;
  }

  void print(pq_container) {
    std::cout << "[NUMA_START]" << std::endl;
    /*for (int i=0; i < numa_pqs.size(); ++i) {
      std::cout << "{DOMAIN:" << i << "}" << std::endl;
      work_item wi;
      while(numa_pqs[i]->pop(wi, 0)) {
        std::cout << "(" <<  std::get<0>(wi) << ", "
                  << std::get<1>(wi) << "), "; 
      }
      std::cout << std::endl;
      std::cout << "{DOMAIN:" << i << "-END}" << std::endl;
      }*/
    std::cout << "[NUMA_END]" << std::endl;    
  }

  void clear(buffer_container) {
    // Note : does not clear the actual numa vector
    for (int i=0; i < numa_buckets.size(); ++i) {
      numa_buckets[i]->clear();
    }
  }

  // TODO revisit clear functions and make sure only append buffers are
  // cleared
  void clear(pq_container) {
    // Note : does not clear the actual numa vector
    for (int i=0; i < numa_pqs.size(); ++i) {
      assert(numa_pqs[i]->empty());
    }
  }

  bool empty(int tid, pq_container) {
    int numanode = this->rt.find_numa_node(tid);
    return numa_pqs[numanode]->empty(tid);
  }

  bool empty(int tid, buffer_container) {
    int numanode = this->rt.find_numa_node(tid);
    return numa_buckets[numanode]->empty(tid);
  }
  
public:
  eagm_node_structure(const work_item& _rep,
                      const eagm_config_t& _config,
                      runtime_t& _rt) : eagm_node_base<work_item,
                                                       current_level_trait>(_rep,
                                                                            _config,
                                                                            _rt){
    numa_container_t c;
    initialize(c);
  }
  
  void push_back(const work_item& wi, int tid) {
    numa_container_t c;
    push_back(wi, tid, c);
  }

  void print() {
    numa_container_t c;
    print(c);
  }

  void clear() {
    numa_container_t c;
    clear(c);
  }

  void process(int tid) {
    numa_container_t c;
    process(tid, c);
  }

  bool empty(int tid) {
    numa_container_t c;
    return empty(tid, c);
  }

private:
  all_numa_buckets_t numa_buckets;
  all_numa_pqs_t numa_pqs;
};

}}}      
#endif
