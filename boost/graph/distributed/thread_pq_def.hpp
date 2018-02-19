// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_GRAPH_THREAD_PRIORITY_Q_DEF
#define BOOST_GRAPH_THREAD_PRIORITY_Q_DEF

#include <queue>

namespace boost {
namespace graph {
namespace distributed {

//==================== The default thread priority queue ==========================
template<typename vertex_distance, typename Compare>
struct thread_priority_queue {

  typedef std::priority_queue<vertex_distance, std::vector<vertex_distance>, Compare> DefaultPriorityQueueType;
  
  thread_priority_queue(const thread_priority_queue& npq) {
  }

  thread_priority_queue(int threads) {
    thread_pqs.resize(threads);
    for(int i=0; i < threads; ++i) {
      thread_pqs[i] = new DefaultPriorityQueueType();
    }
  }

  ~thread_priority_queue() {
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
    if (thread_pqs[tid]->empty())
      return false;

    p = thread_pqs[tid]->top();
    thread_pqs[tid]->pop();
    return true;
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

// Thread priority queue generator
struct thread_priority_queue_gen {
  template<typename vertex_distance, typename Compare>
  struct queue {
    typedef thread_priority_queue<vertex_distance, Compare> type;
  };
};

}
}
}
#endif
