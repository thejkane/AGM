#ifndef PBGL2_SPIN_LOCK
#define PBGL2_SPIN_LOCK

//reference -- http://www.boost.org/doc/libs/1_54_0/doc/html/atomic/usage_examples.html#boost_atomic.usage_examples.example_spinlock

#include <boost/atomic.hpp>

class spinlock {
private:
  typedef enum {Locked, Unlocked} LockState;
  boost::atomic<LockState> state_;

public:
  spinlock() : state_(Unlocked) {}

  void lock()
  {
    while (state_.exchange(Locked, boost::memory_order_acquire) == Locked) {
      /* busy-wait */
    }
  }
  void unlock()
  {
    state_.store(Unlocked, boost::memory_order_release);
  }
};

template<typename Lock>
class vertex_state_lock {
private:
  Lock _lock;

public:
  void lock() {
    _lock.lock();
  }

  void unlock() {
    _lock.unlock();
  }
};

#endif
