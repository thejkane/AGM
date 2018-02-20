// Copyright 2010-2013 The Trustees of Indiana University.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#ifndef AMPLUSPLUS_DETAIL_THREAD_SUPPORT_HPP
#define AMPLUSPLUS_DETAIL_THREAD_SUPPORT_HPP

#include <boost/thread/locks.hpp>
#include <boost/thread/tss.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/assert.hpp>

#ifdef AMPLUSPLUS_SINGLE_THREADED
#include <boost/signals2/dummy_mutex.hpp>
namespace amplusplus {
  namespace detail {
    typedef boost::signals2::dummy_mutex recursive_mutex;
    typedef boost::signals2::dummy_mutex mutex;

    class barrier {
      public:
      barrier(unsigned int count) {(void)count; BOOST_ASSERT (count == 1);}
      bool wait() {return true;}
    };
    inline void do_pause() {}
  }
}
#define AMPLUSPLUS_MULTITHREAD(x) /**/
#else
#include <boost/thread.hpp>
#include <pthread.h>
#if defined(__i386) || defined(__x86_64) && !_CRAYC
#include <xmmintrin.h>
#endif
namespace amplusplus {
  namespace detail {
#if 0
    struct mutex {
      pthread_spinlock_t actual_lock;

      mutex() {pthread_spin_init(&actual_lock, PTHREAD_PROCESS_PRIVATE);}
      ~mutex() {pthread_spin_destroy(&actual_lock);}

      void lock() {pthread_spin_lock(&actual_lock);}
      void unlock() {pthread_spin_unlock(&actual_lock);}
    };
    struct recursive_mutex {
      boost::thread_specific_ptr<int> my_lock_count;
      mutex m;

      recursive_mutex(): my_lock_count(), m() {}

      void lock() {
        if (!my_lock_count.get()) my_lock_count.reset(new int(0));
        int& my_lock_count_data = *my_lock_count;
        if (my_lock_count_data++ == 0) m.lock();
      }

      void unlock() {
        BOOST_ASSERT (my_lock_count.get());
        int& my_lock_count_data = *my_lock_count;
        if (--my_lock_count_data == 0) m.unlock();
      }
    };
#endif
    typedef boost::mutex mutex;
    typedef boost::recursive_mutex recursive_mutex;
    typedef boost::barrier barrier;

#if defined(__i386) || defined(__x86_64) && !_CRAYC
    inline void do_pause() {_mm_pause();}
#else
    inline void do_pause() {}
#endif

  }
}
#define AMPLUSPLUS_MULTITHREAD(x) x
#endif

#ifdef AMPLUSPLUS_SINGLE_THREADED
namespace amplusplus {namespace detail {
template <typename T>
class atomic {
  T value;
  public:
  atomic(T x = T()): value(x) {}
  T load() const {return value;}
  void store(T x) {value = x;}
  T exchange(T x) {T old_value = value; value = x; return old_value;}
  bool compare_exchange_strong(T& old_value, T new_value) {if (value == old_value) {value = new_value; return true;} else {old_value = value; return false;}}
  bool compare_exchange_weak(T& old_value, T new_value) {if (value == old_value) {value = new_value; return true;} else {old_value = value; return false;}}
  T fetch_add(T x) {value += x; return value - x;}
  T fetch_sub(T x) {value -= x; return value + x;}
  T fetch_or(T x) {T old_value = value; value |= x; return old_value;}
  T fetch_and(T x) {T old_value = value; value &= x; return old_value;}
  atomic& operator++() {fetch_add(1); return *this;}
  atomic& operator--() {fetch_add(-1); return *this;}
  atomic& operator+=(T x) {value += x; return *this;}
  atomic& operator-=(T x) {value -= x; return *this;}
};
}}
#elif 1
#ifdef AMPLUSPLUS_BUILTIN_ATOMICS
#include <atomic>
#else
#include <boost/atomic.hpp>
#endif
// #include <cstdatomic>
namespace amplusplus {namespace detail {
#ifdef AMPLUSPLUS_BUILTIN_ATOMICS
using std::atomic;
#else
using boost::atomic;
#endif
}}
#elif 1
namespace amplusplus {namespace detail {

template <typename T, typename = void> struct atomics_supported: boost::mpl::false_ {};

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_1
template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 1>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_2
template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 2>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 4>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_8
template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 8>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_16
template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 16>::type>: boost::mpl::true_ {};
#endif

template <typename T>
class atomic {
  BOOST_MPL_ASSERT((atomics_supported<T>));

  volatile T value;
  public:
  atomic(T x = T()): value(x) {}
  T load() const {__sync_synchronize(); T val = value; __sync_synchronize(); return val;}
  void store(T x) {value = x; __sync_synchronize();}
  T exchange(T x) {
    T prev = 0, prev_old = 0;
    do {
      prev_old = prev;
      prev = __sync_val_compare_and_swap(&value, prev_old, x);
    } while (prev != prev_old);
    return prev;
  }
  bool compare_exchange_strong(T& old_value, T new_value) {
    T prev = __sync_val_compare_and_swap(&value, old_value, new_value);
    if (prev == old_value) return true;
    old_value = prev;
    return false;
  }
  bool compare_exchange_weak(T& old_value, T new_value) {return compare_exchange_strong(old_value, new_value);}
  atomic<T>& operator+=(T x) {fetch_add(x); return *this;}
  atomic<T>& operator-=(T x) {fetch_sub(x); return *this;}
  T fetch_add(T x) {T val = __sync_fetch_and_add(&value, x); return val;}
  T fetch_sub(T x) {T val = __sync_fetch_and_sub(&value, x); return val;}
  T fetch_or(T x) {T val = __sync_fetch_and_or(&value, x); return val;}
  T fetch_and(T x) {T val = __sync_fetch_and_and(&value, x); return val;}
  atomic& operator++() {*this += 1; return *this;}
  atomic& operator--() {*this -= 1; return *this;}
};
}}
#else
// Lock-based version for thread debugging tools
namespace amplusplus {namespace detail {
template <typename T>
class atomic: boost::noncopyable {
  mutable boost::mutex lock;
  T value;
  public:
  atomic(T x = T()): lock(), value(x) {}
  T load() const {boost::lock_guard<boost::mutex> l(lock); return value;}
  void store(T x) {boost::lock_guard<boost::mutex> l(lock); value = x;}
  T exchange(T x) {
    boost::lock_guard<boost::mutex> l(lock);
    std::swap(value, x);
    return x;
  }
  bool compare_exchange_strong(T& old_value, T new_value) {
    boost::lock_guard<boost::mutex> l(lock);
    if (value == old_value) {
      value = new_value;
      return true;
    } else {
      old_value = value;
      return false;
    }
  }
  bool compare_exchange_weak(T& old_value, T new_value) {return compare_exchange_strong(old_value, new_value);}
  T fetch_add(T x) {boost::lock_guard<boost::mutex> l(lock); T old_value = value; value += x; return old_value;}
  T fetch_sub(T x) {return fetch_add(-x);}
  T fetch_or(T x) {boost::lock_guard<boost::mutex> l(lock); T old_value = value; value |= x; return old_value;}
  T fetch_and(T x) {boost::lock_guard<boost::mutex> l(lock); T old_value = value; value &= x; return old_value;}
  atomic& operator++() {fetch_add(1); return *this;}
  atomic& operator--() {fetch_add(-1); return *this;}
  atomic& operator+=(T x) {fetch_add(x); return *this;}
  atomic& operator-=(T x) {fetch_sub(x); return *this;}
};
}}
#endif

namespace amplusplus {
  namespace detail {
    extern __thread int internal_thread_id;

    static inline int get_thread_id() {
      BOOST_ASSERT (internal_thread_id != -1); // Ensure that it has been set
      return internal_thread_id;
    }

    class push_thread_id_obj {
      int old_id;

      public:
      push_thread_id_obj(int new_id) // Can't be explicit
        : old_id(internal_thread_id)
        {internal_thread_id = new_id;}
      ~push_thread_id_obj() {internal_thread_id = old_id;}
      operator bool() const {return false;} // For use in if statements
    };

    // NOTE : This will break if return value optimizations (RVO) are disabled. When RVO is disabled an additional
    // object will get created inside the if condition and will be deleted, leaving internal_thread_id at -1.
    // Therefore, make sure code is compiled RVO enabled (default) when using this macro.
    #define AMPLUSPLUS_WITH_THREAD_ID(id) if (::amplusplus::detail::push_thread_id_obj BOOST_PP_CAT(tidobj_, __LINE__) = (id)) {} else
  }
}

#endif // AMPLUSPLUS_DETAIL_THREAD_SUPPORT_HPP
