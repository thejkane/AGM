#ifndef __AGM_PRIORITY_Q_H__
#define __AGM_PRIORITY_Q_H__

// LibCDS stuff
//#include <cds/init.h>       // for cds::Initialize and cds::Terminate
//#include <cds/gc/hp.h>      // for cds::HP (Hazard Pointer) SMR
#include <cds/container/fcpriority_queue.h>

namespace boost { namespace graph { namespace agm {
      
template<typename work_item, typename Compare>
struct concurrent_priority_queue {

  typedef std::priority_queue<work_item, std::vector<work_item>, Compare> DefaultPriorityQueueType;
  typedef cds::container::FCPriorityQueue<work_item, DefaultPriorityQueueType> NodePriorityQueueType;

  concurrent_priority_queue(const concurrent_priority_queue& npq) {
  }

  // capacity ignored
  concurrent_priority_queue(int threads, unsigned long cap=0) {
    // Initialize Hazard Pointer singleton
    // cds::gc::HP hpGC;

    // Attach thread
    //    cds::threading::Manager::attachThread();
  }

  void put(const work_item& p, int tid) {
    node_pq.push(p);
  }

  bool pop(work_item& p, int tid) {
    return node_pq.pop(p);
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

}}}      
#endif
