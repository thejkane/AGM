// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_GRAPH_AGM_PRIORITY_Q_DEFS
#define BOOST_GRAPH_AGM_PRIORITY_Q_DEFS


#include <queue>
#include <algorithm> // for std::min, std::max
#include <iostream>
#include <atomic>
#include <tuple>
#include <climits>

// LibCDS stuff
#include <cds/init.h>       // for cds::Initialize and cds::Terminate
#include <cds/gc/hp.h>      // for cds::HP (Hazard Pointer) SMR
#include <cds/container/fcpriority_queue.h>
#include <cds/container/mspriority_queue.h>
//#include "ellen_bintree_pqueue.h"
#include "thread_pq_def.hpp"

// NUMA32
#ifndef __APPLE__
#include <numaif.h>
#include <numa.h>
#include <unistd.h>

#endif

namespace boost {
namespace graph {
namespace distributed {

// The multi priority queue
template<typename vertex_distance, typename Compare>
struct multi_std_priority_queue {

#define PQ_DEFINED 2

  typedef multi_std_priority_queue<vertex_distance, Compare> Self;
  typedef std::priority_queue<vertex_distance, std::vector<vertex_distance>, Compare> DefaultPriorityQueueType;

  int put_index = 0;
  int pop_index = 0;
  DefaultPriorityQueueType priority_q[PQ_DEFINED];

  multi_std_priority_queue() {}

  void push(const vertex_distance& p) {
    priority_q[put_index].push(p);
    put_index = (++put_index)%PQ_DEFINED;
  }

  bool pop(vertex_distance& topv) {
    if (!priority_q[pop_index].empty()) {
      topv = priority_q[pop_index].top();
      priority_q[pop_index].pop();
      pop_index = (++pop_index)%PQ_DEFINED;
      return true;
    }

    for (int i=1; i < PQ_DEFINED; ++i) {
      int index = (i+pop_index)%PQ_DEFINED;
      if (!priority_q[index].empty()) {
	pop_index = index;
	topv = priority_q[pop_index].top();
	priority_q[pop_index].pop();
	return true;
      }
    }

    return false;

  }

  bool empty() const {
    if (!priority_q[pop_index].empty()) {
      return false;
    }

    for (int i=1; i < PQ_DEFINED; ++i) {
      int index = (i+pop_index)%PQ_DEFINED;
      if (!priority_q[index].empty()) {
	const_cast<Self*>(this)->pop_index = index;
	return false;
      }
    }

    return true;
  }

  size_t size() const {

    size_t size=0;
    for(int i=0; i<PQ_DEFINED; ++i) {
      size += priority_q[i].size();
    }
    
    return size;
  }

};

//==================== The default thread priority queue ==========================
template<typename vertex_distance, typename Compare>
struct multi_priority_queue {

  typedef multi_std_priority_queue<vertex_distance, Compare> DefaultPriorityQueueType;
  

  multi_priority_queue(int threads) {
    thread_pqs.resize(threads);
    for(int i=0; i < threads; ++i) {
      thread_pqs[i] = new DefaultPriorityQueueType();
    }
  }

  ~multi_priority_queue() {
    typename std::vector<DefaultPriorityQueueType*>::iterator 
      ite = thread_pqs.begin();
    for(; ite != thread_pqs.end(); ++ite) {
      delete (*ite);
    }
  }

  void initialize(int _tid, bool ignore) {}

  void allocate_pqs() {}

  void put(const vertex_distance& p, int tid) {
    thread_pqs[tid]->push(p);
  }

  bool pop(vertex_distance& p, int tid) {
    return thread_pqs[tid]->pop(p);
  }

  vertex_distance top(int tid) {
    assert(false);
  }

  size_t size(int tid) const {
    return thread_pqs[tid]->size();
  }

  size_t size() const {
    //TODO 
    assert(false);
  }

  void clear(int tid) {
    assert(thread_pqs[tid]->empty());
  }

  void clear() {
    // TODO
    assert(false);
  }

  bool empty(int tid) {
    return thread_pqs[tid]->empty();
  }

  bool empty() {
    // this call is not thread safe
    bool emp = true;
    for (int i=0; i<thread_pqs.size(); ++i) {
      emp = emp && empty(i);
    }
    return emp;
  }

private:
  std::vector<DefaultPriorityQueueType*> thread_pqs;
};


// Multi priority queue generator
struct multi_priority_queue_gen {
  template<typename vertex_distance, typename Compare>
  struct queue {
    typedef multi_priority_queue<vertex_distance, Compare> type;
  };
};

//================= The default node priority queue ================================
template<typename vertex_distance, typename Compare>
struct node_priority_queue {

  typedef std::priority_queue<vertex_distance, std::vector<vertex_distance>, Compare> DefaultPriorityQueueType;
  typedef cds::container::FCPriorityQueue<vertex_distance, DefaultPriorityQueueType> NodePriorityQueueType;

  node_priority_queue(const node_priority_queue& npq) {
  }

  // capacity ignored
  node_priority_queue(int threads, unsigned long cap=0) {
    // Initialize Hazard Pointer singleton
    // cds::gc::HP hpGC;

    // Attach thread
    //    cds::threading::Manager::attachThread();
  }

  void initialize(int _tid, bool ignore) {}

  void allocate_pqs() {}

  void put(const vertex_distance& p, int tid) {
    node_pq.push(p);
  }

  bool pop(vertex_distance& p, int tid) {
    return node_pq.pop(p);
  }

  vertex_distance top(int tid) {
    assert(false);
  }

  size_t size(int tid) const {
    return node_pq.size();
  }

  size_t size() const {
    return size(0); // tid doesnt matter
  }

  void clear() {
    node_pq.clear();
  }

  bool empty(int tid) {
    return node_pq.empty();
  }

  bool empty() {
    return empty(0);
  }


private:
  NodePriorityQueueType node_pq;
};

// Node priority queue generator
struct node_priority_queue_gen {
  template<typename vertex_distance, typename Compare>
  struct queue {
    typedef node_priority_queue<vertex_distance, Compare> type;
  };
};



template<typename vertex_distance_data, typename Compare>
struct reverse_compare {
private:
  Compare c;
public:
  inline bool operator()(const vertex_distance_data& vd1, const vertex_distance_data& vd2) {
    return c(vd2, vd1);
  }
};


#ifndef RAVEN
#define DEFAULT_HEAP_SIZE (1UL << 29)
#warning "========================== Heap Size 1 << 29 ========================="
#else
#warning "========================== Heap Size 32 << 16 ========================="
#define DEFAULT_HEAP_SIZE (32 << 16)
#endif

  //================== MS PQ Based Node Q =======================================
template<typename vertex_distance, typename Compare>
struct ms_node_priority_queue {
  // we need this bcos MSPriorityQueue is ordering in reverse
  typedef reverse_compare<vertex_distance, Compare> ReverseCompare;

  struct traits_MSPriorityQueue_dyn_cmp : public
  cds::container::mspriority_queue::make_traits <
    cds::opt::buffer< cds::opt::v::dynamic_buffer< char > >
    , cds::opt::compare < ReverseCompare >, cds::opt::less < std::greater<vertex_distance> >
    > ::type
  {};

  // for scale 20 - 16
  // for scale 15 - 10 (barely)

  typedef cds::container::MSPriorityQueue<vertex_distance, 
					  traits_MSPriorityQueue_dyn_cmp> NodePriorityQueueType;

  ms_node_priority_queue(const ms_node_priority_queue& npq) : node_pq(npq.capacity) {
  }

  ms_node_priority_queue(int threads, unsigned long cap=DEFAULT_HEAP_SIZE) : capacity(cap), 
							   node_pq(cap) {
    
    // Initialize Hazard Pointer singleton
    // cds::gc::HP hpGC;

    // Attach thread
    //    cds::threading::Manager::attachThread();
  }

  void initialize(int _tid, bool ignore) {}

  void allocate_pqs() {}

  void put(const vertex_distance& p, int tid) {
    if (!node_pq.push(p)) {
      std::cerr << "Not enough space ms_node_priority_queue !, Total queue size : " << node_pq.size() << std::endl;
      assert(false);
      exit(-1);
    }
  }

  bool pop(vertex_distance& p, int tid) {
    return node_pq.pop(p);
  }

  vertex_distance top(int tid) {
    assert(false);
  }

  size_t size(int tid) const {
    return node_pq.size();
  }

  size_t size() const {
    return size(0); // tid doesnt matter
  }

  void clear() {
    node_pq.clear();
  }

  bool empty(int tid) {
    return node_pq.empty();
  }

  bool empty() {
    return empty(0);
  }


private:
  NodePriorityQueueType node_pq;
  unsigned long capacity;
};

// Node priority queue generator
struct ms_node_priority_queue_gen {
  template<typename vertex_distance, typename Compare>
  struct queue {
    typedef ms_node_priority_queue<vertex_distance, Compare> type;
  };
};



  //================== Ellenbin Based Node Q =======================================
  /*template<typename vertex_distance, typename Compare>
struct ellen_bin_priority_queue {
  // we need this bcos MSPriorityQueue is ordering in reverse
  typedef reverse_compare<vertex_distance, Compare> ReverseCompare;
  typedef vertex_distance KeyType;
  
  struct dist_key_extractor {
    void operator ()( KeyType& dest, vertex_distance const& src ) {
      dest = src; // key is distance
    }
  };


  struct key_comp {
    bool operator()( vertex_distance const& v1, vertex_distance const& v2 ) const { 
      std::cout << "v1.f1 : " << v1.first << " v1.f2 :" << v1.second
		<< " v2.f1 : " << v2.first << " v2.f2 :" << v2.second << std::endl;
      if (v1.second < v2.second)
	return false;
      else if (v1.second == v2.second) {
	return v1.first < v2.first;
      } else
	return true;
    }


    // Support comparing KeyType and char const 

  };

 
  struct traits_EllenBinTree_min :
    public cds::container::ellen_bintree::make_set_traits<
    cds::container::ellen_bintree::key_extractor< dist_key_extractor >
    , cds::opt::compare < key_comp >
    >::type
  {};

  typedef pqueue::EllenBinTreePQueue< cds::gc::HP, 
			      KeyType, 
			      vertex_distance, 
			      traits_EllenBinTree_min, false > NodePriorityQueueType;


  ellen_bin_priority_queue(const ellen_bin_priority_queue& npq) {
  }

  ellen_bin_priority_queue(int threads) {
    
    // Initialize Hazard Pointer singleton
    // cds::gc::HP hpGC;

    // Attach thread
    //    cds::threading::Manager::attachThread();
  }

  void initialize(int _tid, bool ignore) {}

  void allocate_pqs() {}

  void put(const vertex_distance& p, int tid) {
    node_pq.push(p);
  }

  bool pop(vertex_distance& p, int tid) {
    return node_pq.pop(p);
  }

  vertex_distance top(int tid) {
    assert(false);
  }

  size_t size(int tid) const {
    return node_pq.size();
  }

  size_t size() const {
    return size(0); // tid doesnt matter
  }

  void clear() {
    node_pq.clear();
  }

  bool empty(int tid) {
    return node_pq.empty();
  }

  bool empty() {
    return empty(0);
  }


private:
  NodePriorityQueueType node_pq;
};

// Node priority queue generator
struct ellen_bin_priority_queue_gen {
  template<typename vertex_distance, typename Compare>
  struct queue {
    typedef ellen_bin_priority_queue<vertex_distance, Compare> type;
  };
};
  */

boost::mutex created_mutex;

int pin(int core) {
#ifndef __APPLE__
  int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
  if (core < 0 || core >= num_cores) {
    std::cerr << "Invalid core " << core << std::endl;
    return 1;
  }

  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core, &cpuset);
  // In Linux underneath boost is using pthread. May not compatible
  // with other platforms
  pthread_t current_thread = pthread_self();    
  return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
#else
  return 0;
#endif
}

//==================== The default NUMA priority queue ==========================
template<typename vertex_distance, typename Compare>
struct numa_priority_queue {

  typedef std::priority_queue<vertex_distance, std::vector<vertex_distance>, Compare> DefaultPriorityQueueType;
  typedef cds::container::FCPriorityQueue<vertex_distance, DefaultPriorityQueueType> NodePriorityQueueType;
  typedef std::vector<NodePriorityQueueType*> PriorityQueueVector_t;

  numa_priority_queue(const numa_priority_queue& npq) {
    std::cout << "Calling copy constructor ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    assert(false);
  }

  numa_priority_queue(int nthreads, unsigned long cap=0) {
#ifndef __APPLE__
    // check whether NUMA available
    // TODO may need to go to the constructor
    if (numa_available() < 0) {
      std::cerr << "System does not support numa\n" << std::endl;
      exit(1);
    }
#ifdef PRINT_DEBUG
    std::cout << "Instantiating numa queue for threads : " << nthreads 
	      << std::endl;
#endif
    thread_node_map = new int[nthreads];
    int max = numa_max_node();
    pq_vector.resize(max+1);
    //    created = new std::atomic<int>[max+1];
    created = new int[max+1];
    for(int i=0; i < (max+1); ++i) {
      //created[i] = ATOMIC_VAR_INIT(0);
      created[i] = 0;
    }
#endif    
  }

  ~numa_priority_queue() {
    delete[] thread_node_map;
    delete[] created;
    
#ifdef PRINT_DEBUG
    std::set<int> deleted_nodes;
    std::set<int>::iterator it;
#endif
    for (int i=0; i < used_nodes.size(); ++i) {
#ifdef PRINT_DEBUG
      it = deleted_nodes.find(used_nodes[i]);
      assert(it == deleted_nodes.end() && "Same node cannot be in the set again !!");
#endif

      delete pq_vector[(used_nodes[i])];

#ifdef PRINT_DEBUG
      deleted_nodes.insert(used_nodes[i]);
#endif
    }
  
    pq_vector.clear();
    used_nodes.clear();
  }

  void allocate_pqs() {
    for (int i=0; i < used_nodes.size(); ++i) {
      //      pq_vector[(used_nodes[i])] = new NodePriorityQueueType();
    }
  }

  // call per each thread
  void initialize(int _tid, bool pin_thread = true) {

    if (pin_thread) {
      // first pin the thread
      if (pin(_tid)) {
	std::cerr << "Error pining thread to cpu ...exiting..." << std::endl;
	exit(1);
      }
    }

    // find the numa node of this thread
    // => find the cpu of this thread
    // => find the numa node
#ifndef __APPLE__    
    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    assert(cpu == _tid);
    thread_node_map[_tid] = node;
#endif    

    // Set policy saying that current thread 
    // can only allocate memory from the numa node
    // that core belongs to current thread
    // we may not need this, cos default policy is to
    // allocate in local
    /*unsigned long mask;
      mask = 1L << node;
      // TODO verify with Marcin
      if (set_mempolicy(MPOL_BIND, &mask, sizeof(mask)*8)) {
      std::cout << "Error setting memory policy " << std::endl;
      exit(1);
      }*/

    //    if (!std::atomic_fetch_add(&created[node], 1)) {
    {
      boost::unique_lock<boost::mutex> scoped_lock(created_mutex);
      if (created[node] == 0) {
#ifdef PRINT_DEBUG
	auto result = std::find(std::begin(used_nodes), std::end(used_nodes), node);
	assert (result == std::end(used_nodes));
#endif
	pq_vector[node] = new NodePriorityQueueType();
	used_nodes.push_back(node);      
	created[node] = 1;
      }
    }
  }

  void put(const vertex_distance& p, int tid) {
    NodePriorityQueueType* node_pq = find_queue(tid);
    node_pq->push(p);
  }

  bool pop(vertex_distance& p, int tid) {
    NodePriorityQueueType* node_pq = find_queue(tid);
    bool ret = node_pq->pop(p);
    return ret;
  }

  bool empty(int tid) const {
    NodePriorityQueueType* node_pq = find_queue(tid);
    bool empty = node_pq->empty();
    return empty;
  }

  vertex_distance top(int tid) {
    assert(false);
  }

  size_t size(int tid) const {
    NodePriorityQueueType* node_pq = find_queue(tid);
    return node_pq->size();
  }

  void clear(int tid) const {
    NodePriorityQueueType* node_pq = find_queue(tid);
    node_pq->clear();
  }

  bool empty() {
    // this call is not thread safe
    bool emp = true;
    for (int i=0; i < used_nodes.size(); ++i) {
      if (!pq_vector[(used_nodes[i])]->empty()) {
	emp = false;
	break;
      }
    }
    return emp;
  }


private:
  int* thread_node_map;
  PriorityQueueVector_t pq_vector;
  //std::atomic<int>* created;
  int* created;
  // used nodes change depending upon thread configuration
  std::vector<int> used_nodes; 

  NodePriorityQueueType* find_queue(int _tid) const {
    // find the node id for thread id
    int nodeid = thread_node_map[_tid];
#ifdef PRINT_DEBUG
    std::cout << "Finding q for tid : " << _tid << " and node id : " << nodeid
	      << std::endl;
#endif
    // find the numa node queue
    NodePriorityQueueType* q = pq_vector[nodeid];
    assert(q != NULL);
#ifdef PRINT_DEBUG
    std::cout << "Pushing to node id : " << nodeid << std::endl;
#endif
    return q;
  }

};

// Node priority queue generator
struct numa_priority_queue_gen {
  template<typename vertex_distance_t, typename Compare>
  struct queue {
    typedef numa_priority_queue<vertex_distance_t, Compare> type;
  };
};


//==================== The default NUMA priority queue ==========================
template<typename vertex_distance, typename Compare>
struct ms_numa_priority_queue {

  // we need this bcos MSPriorityQueue is ordering in reverse
  typedef reverse_compare<vertex_distance, Compare> ReverseCompare;

  struct traits_MSPriorityQueue_dyn_cmp : public
  cds::container::mspriority_queue::make_traits <
    cds::opt::buffer< cds::opt::v::dynamic_buffer< char > >
    , cds::opt::compare < ReverseCompare >, cds::opt::less < std::greater<vertex_distance> >
    > ::type
  {};


  typedef cds::container::MSPriorityQueue<vertex_distance, 
					  traits_MSPriorityQueue_dyn_cmp> NodePriorityQueueType;

  typedef std::vector<NodePriorityQueueType*> PriorityQueueVector_t;

  ms_numa_priority_queue(const ms_numa_priority_queue& npq) {
    std::cout << "Calling copy constructor ~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    assert(false);
  }

  ms_numa_priority_queue(int nthreads, unsigned long cap=DEFAULT_HEAP_SIZE) : capacity(cap) {
#ifndef __APPLE__    
    // check whether NUMA available
    // TODO may need to go to the constructor
    if (numa_available() < 0) {
      std::cerr << "System does not support numa\n" << std::endl;
      exit(1);
    }

#ifdef PRINT_DEBUG
    std::cout << "Instantiating numa queue for threads : " << nthreads 
	      << std::endl;
#endif
    thread_node_map = new int[nthreads];
    int max = numa_max_node();
    pq_vector.resize(max+1);
    //    created = new std::atomic<int>[max+1];
    created = new int[max+1];
    for(int i=0; i < (max+1); ++i) {
      //created[i] = ATOMIC_VAR_INIT(0);
      created[i] = 0;
    }
#endif    
  }

  ~ms_numa_priority_queue() {
    delete[] thread_node_map;
    delete[] created;
    
#ifdef PRINT_DEBUG
    std::set<int> deleted_nodes;
    std::set<int>::iterator it;
#endif
    for (int i=0; i < used_nodes.size(); ++i) {
#ifdef PRINT_DEBUG
      it = deleted_nodes.find(used_nodes[i]);
      assert(it == deleted_nodes.end() && "Same node cannot be in the set again !!");
#endif

      delete pq_vector[(used_nodes[i])];

#ifdef PRINT_DEBUG
      deleted_nodes.insert(used_nodes[i]);
#endif
    }
  
    pq_vector.clear();
    used_nodes.clear();
  }

  void allocate_pqs() {
    assert(false);
  }

  // call per each thread
  void initialize(int _tid, bool pin_thread = true) {

    if (pin_thread) {
      // first pin the thread
      if (pin(_tid)) {
	std::cerr << "Error pining thread to cpu ...exiting..." << std::endl;
	exit(1);
      }
    }

    // find the numa node of this thread
    // => find the cpu of this thread
    // => find the numa node
#ifndef __APPLE__    
    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    assert(cpu == _tid);
    thread_node_map[_tid] = node;
#endif    

    // Set policy saying that current thread 
    // can only allocate memory from the numa node
    // that core belongs to current thread
    // we may not need this, cos default policy is to
    // allocate in local
    /*unsigned long mask;
      mask = 1L << node;
      // TODO verify with Marcin
      if (set_mempolicy(MPOL_BIND, &mask, sizeof(mask)*8)) {
      std::cout << "Error setting memory policy " << std::endl;
      exit(1);
      }*/

    //    if (!std::atomic_fetch_add(&created[node], 1)) {
    {
      boost::unique_lock<boost::mutex> scoped_lock(created_mutex);
      if (created[node] == 0) {
#ifdef PRINT_DEBUG
	auto result = std::find(std::begin(used_nodes), std::end(used_nodes), node);
	assert (result == std::end(used_nodes));
#endif
	pq_vector[node] = new NodePriorityQueueType(capacity/pq_vector.size());
	used_nodes.push_back(node);      
	created[node] = 1;
      }
    }
  }

  void put(const vertex_distance& p, int tid) {
    NodePriorityQueueType* node_pq = find_queue(tid);
    node_pq->push(p);
  }

  bool pop(vertex_distance& p, int tid) {
    NodePriorityQueueType* node_pq = find_queue(tid);
    bool ret = node_pq->pop(p);
    return ret;
  }

  bool empty(int tid) const {
    NodePriorityQueueType* node_pq = find_queue(tid);
    bool empty = node_pq->empty();
    return empty;
  }

  vertex_distance top(int tid) {
    assert(false);
  }

  size_t size(int tid) const {
    NodePriorityQueueType* node_pq = find_queue(tid);
    return node_pq->size();
  }

  void clear(int tid) const {
    NodePriorityQueueType* node_pq = find_queue(tid);
    node_pq->clear();
  }

  bool empty() {
    // this call is not thread safe
    bool emp = true;
    for (int i=0; i < used_nodes.size(); ++i) {
      if (!pq_vector[(used_nodes[i])]->empty()) {
	emp = false;
	break;
      }
    }
    return emp;
  }


private:
  int* thread_node_map;
  PriorityQueueVector_t pq_vector;
  //std::atomic<int>* created;
  int* created;
  // used nodes change depending upon thread configuration
  std::vector<int> used_nodes; 
  // capacity
  unsigned long capacity;

  NodePriorityQueueType* find_queue(int _tid) const {
    // find the node id for thread id
    int nodeid = thread_node_map[_tid];
#ifdef PRINT_DEBUG
    std::cout << "Finding q for tid : " << _tid << " and node id : " << nodeid
	      << std::endl;
#endif
    // find the numa node queue
    NodePriorityQueueType* q = pq_vector[nodeid];
    assert(q != NULL);
#ifdef PRINT_DEBUG
    std::cout << "Pushing to node id : " << nodeid << std::endl;
#endif
    return q;
  }

};

// Node priority queue generator
struct ms_numa_priority_queue_gen {
  template<typename vertex_distance_t, typename Compare>
  struct queue {
    typedef ms_numa_priority_queue<vertex_distance_t, Compare> type;
  };
};



} // end namespace distributed 
} // end namespace graph
} // end namespace boost

#endif
