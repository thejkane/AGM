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

#ifndef AMPLUSPLUS_MESSAGE_QUEUE_HPP
#define AMPLUSPLUS_MESSAGE_QUEUE_HPP

#include <am++/traits.hpp>
#include <am++/detail/thread_support.hpp>
#include <boost/config.hpp>
#include <boost/bind.hpp>
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/assert.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/intrusive/slist.hpp>
#include <list>
#include <stdio.h>

namespace amplusplus {

class scheduler;

// #define AMPLUSPLUS_USE_STD_LIST // For easier debugging

#ifdef AMPLUSPLUS_USE_STD_LIST

struct task_hook_type {};

#else // Use Boost.Intrusive

struct task_tag {};
typedef boost::intrusive::slist_base_hook<boost::intrusive::tag<task_tag> > task_hook_type;

#endif

struct task_base: task_hook_type {
  virtual int /* task_result */ operator()(scheduler&) const = 0;
  virtual ~task_base() {}
};

template <typename F>
class task_impl: public task_base {
  F f;
  public:
  explicit task_impl(F f): f(AMPLUSPLUS_MOVE(f)) {}
  virtual int /* task_result */ operator()(scheduler& s) const {return f(s);}
};

typedef task_base* task;

class scheduler: boost::noncopyable {
  public:
  enum task_result {tr_remove_from_queue, tr_idle, tr_busy, tr_busy_and_finished};
  typedef task function_type;

  private:
  amplusplus::detail::mutex lock;
#ifdef AMPLUSPLUS_USE_STD_LIST
  typedef std::list<task_base*> run_queue_type;
#else
  typedef boost::intrusive::slist<
            task_base,
            boost::intrusive::base_hook<task_hook_type>,
            boost::intrusive::constant_time_size<false>,
            boost::intrusive::cache_last<true>
          > run_queue_type;
#endif
  run_queue_type run_queue;
  boost::thread_specific_ptr<unsigned int> reentry_count; // Only counts reentries that should disable running handlers.
  // run_queue_type idle_tasks;

#ifndef BOOST_NO_DELETED_FUNCTIONS
  scheduler(const scheduler&) = delete;
  scheduler(scheduler&&) = delete;
#endif

  struct delete_task {void operator()(task t) const {delete t;}};

  public:
  scheduler() {}
  ~scheduler() {
    this->run_until(boost::bind(&run_queue_type::empty, &run_queue));
    BOOST_ASSERT (run_queue.empty());
    // idle_tasks.clear_and_dispose(delete_task());
  }

  template <typename F, int priority = 0>
  void add_runnable(F f) {
    boost::lock_guard<amplusplus::detail::mutex> l(lock);
#ifdef AMPLUSPLUS_USE_STD_LIST
    if(priority == 0)
      run_queue.push_back(new task_impl<F>(AMPLUSPLUS_MOVE(f)));
    else {
      // std::cout << "!!!!!!!!!!!!!!!!! Sending Priority !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      run_queue.push_front(new task_impl<F>(AMPLUSPLUS_MOVE(f)));
    }
#else
    if(priority == 0)
      run_queue.push_back(*new task_impl<F>(AMPLUSPLUS_MOVE(f)));
    else {
      // std::cout << "!!!!!!!!!!!!!!!!! Sending Priority !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      run_queue.push_front(*new task_impl<F>(AMPLUSPLUS_MOVE(f)));
    }
#endif
    // fprintf(stderr, "Pushing task in %p\n", (void*)pthread_self());
  }

  template <typename F>
  void add_idle_task(F f) {
    boost::lock_guard<amplusplus::detail::mutex> l(lock);
    // idle_tasks.push_back(*new task_impl<F>(f));
    task t = new task_impl<F>(AMPLUSPLUS_MOVE(f));
#ifdef AMPLUSPLUS_USE_STD_LIST
    run_queue.push_back(t);
#else
    run_queue.push_back(*t);
#endif
    // fprintf(stderr, "Pushing idle task %p in %p\n", t, (void*)pthread_self());
  }

  bool run_one_task(run_queue_type& q) {
    bool any_busy = false;
    task t = 0;
    {
      boost::lock_guard<amplusplus::detail::mutex> l(lock);
      if (q.empty()) {
        // fprintf(stderr, "No task found for %p\n", (void*)pthread_self());
        return false;
      }
#ifdef AMPLUSPLUS_USE_STD_LIST
      t = q.front();
#else
      t = &q.front();
#endif
      q.pop_front();
    }
    assert (t != 0);
    // fprintf(stderr, "Running task %p (%s) in %p\n", t, typeid(*t).name(), (void*)pthread_self());
    task_result r = task_result((*t)(*this));
    // fprintf(stderr, "Task result = %d\n", (int)r);
    switch (r) {
      case tr_busy_and_finished: any_busy = true; // Fall through
      case tr_remove_from_queue: delete t; break;
      case tr_busy: any_busy = true; // Fall through
      case tr_idle:
#ifdef AMPLUSPLUS_USE_STD_LIST
        {boost::lock_guard<amplusplus::detail::mutex> l(lock); q.push_back(t);}
#else
        {boost::lock_guard<amplusplus::detail::mutex> l(lock); q.push_back(*t);}
#endif
        break;
      default: abort();
    }

    return any_busy;
  }

  void run_one() {run_one_task(run_queue);}

  template <typename Pred>
  void run_until(Pred pred) {
    while (!pred()) run_one_task(run_queue);
  }

  template <typename Pred>
  void run_until_for_flow_control(Pred pred) {
    if (!reentry_count.get()) reentry_count.reset(new unsigned int(0));
    ++*reentry_count;
    while (!pred()) run_one_task(run_queue);
    --*reentry_count;
  }

  bool should_run_handlers() const {
    return (!reentry_count.get() || *reentry_count == 0);
  }

};

namespace detail {
  template <typename F>
  class delay_t {
    F f;
    scheduler& sched;
    template <typename Binder>
    struct body_t {
      Binder b;
      explicit body_t(const Binder& b): b(b) {}
      scheduler::task_result operator()(scheduler&) const {b(); return scheduler::tr_busy_and_finished;}
    };
    template <typename Binder> static body_t<Binder> body(const Binder& b) {return body_t<Binder>(b);}
    public:
    explicit delay_t(const F& f, scheduler& sched): f(f), sched(sched) {}
    template <typename Val>
    void operator()(const Val& val) const {sched.add_runnable(body(boost::bind(f, val)));}
  };
}

template <typename F>
detail::delay_t<F> delay(const F& f, scheduler& sched) {
  return detail::delay_t<F>(f, sched);
}

template <typename Val>
class message_queue: boost::noncopyable {
  amplusplus::detail::mutex lock;
  std::list<Val> messages;
  std::list<boost::function<void(Val)> > receivers_waiting;
  bool receive_all_active;
  boost::shared_ptr<bool> alive;

  public:
  message_queue(scheduler&): lock(), messages(), receivers_waiting(), receive_all_active(false), alive(new bool(true)) {}
  ~message_queue() {
    // fprintf(stderr, "message_queue::~message_queue %p\n", this);
    BOOST_ASSERT (messages.empty());
    *alive = false;
    // Might have receivers waiting because of resubmits from things like queue copies
  }

#ifndef BOOST_NO_DELETED_FUNCTIONS
  message_queue(const message_queue&) = delete;
#endif

  void send(Val msg) {
    // fprintf(stderr, "%p getting send in %p\n", this, (void*)pthread_self());
    boost::function<void(Val)> k;
    {
      boost::lock_guard<amplusplus::detail::mutex> l(lock);
      if (!receivers_waiting.empty()) {
        BOOST_ASSERT (messages.empty());
        // fprintf(stderr, "%p have receiver, pushing task in %p\n", this, (void*)pthread_self());
        if (receive_all_active && receivers_waiting.size() == 1) {
          k = receivers_waiting.front();
        } else {
          k.swap(receivers_waiting.front());
          receivers_waiting.pop_front();
        }
        goto outside_lock;
      } else {
        // fprintf(stderr, "%p waiting for receiver in %p\n", this, (void*)pthread_self());
        messages.push_back(AMPLUSPLUS_MOVE(msg));
        return;
      }
    }
    outside_lock:
    k(AMPLUSPLUS_MOVE(msg));
  }

  template <typename Iter>
  void send_range(Iter b, Iter e) {
    // fprintf(stderr, "%p getting send in %p\n", this, (void*)pthread_self());
    while (b != e) {
      boost::function<void(Val)> k;
      bool process_rest = false;
      {
        boost::lock_guard<amplusplus::detail::mutex> l(lock);
        if (!receivers_waiting.empty()) {
          BOOST_ASSERT (messages.empty());
          // fprintf(stderr, "%p have receiver, pushing task in %p\n", this, (void*)pthread_self());
          if (receive_all_active && receivers_waiting.size() == 1) {
            k = receivers_waiting.front();
            process_rest = true;
          } else {
            k.swap(receivers_waiting.front());
            receivers_waiting.pop_front();
          }
          goto outside_lock;
        } else {
          // fprintf(stderr, "%p waiting for receiver in %p\n", this, (void*)pthread_self());
          messages.insert(messages.end(), b, e);
          b = e;
          goto end_of_loop;
        }
      }
      outside_lock:
      if (process_rest) {
        for (; b != e; ++b) k(*b);
      } else {
        k(*b++);
      }

      end_of_loop:;
    }
  }

  template <typename K>
  void receive(K k) {
    BOOST_ASSERT (!receive_all_active);
    // fprintf(stderr, "%p getting receive in %p\n", this, (void*)pthread_self());
    boost::optional<Val> msg;
    {
      boost::lock_guard<amplusplus::detail::mutex> l(lock);
      if (!messages.empty()) {
        // fprintf(stderr, "%p have message, pushing task in %p\n", this, (void*)pthread_self());
        BOOST_ASSERT (receivers_waiting.empty());
        msg = AMPLUSPLUS_MOVE(messages.front());
        messages.pop_front();
        goto outside_lock;
      } else {
        // fprintf(stderr, "%p waiting for message in %p\n", this, (void*)pthread_self());
        receivers_waiting.push_back(k);
        return;
      }
    }
    outside_lock:
    AMPLUSPLUS_MOVE(k)(msg.get());
  }

  template <typename K>
  void receive_all(K k) {
    BOOST_ASSERT (!receive_all_active);
    std::list<Val> messages_to_process;
    {
      boost::lock_guard<amplusplus::detail::mutex> l(lock);
      messages_to_process.swap(messages);
    }
    for (typename std::list<Val>::const_iterator i = messages_to_process.begin();
         i != messages_to_process.end(); ++i) {
      k(*i);
    }
    {
      boost::lock_guard<amplusplus::detail::mutex> l(lock);
      receivers_waiting.push_back(AMPLUSPLUS_MOVE(k));
      receive_all_active = true;
    }
  }

  boost::shared_ptr<bool> get_alive() const {return alive;}
};

template <typename Val>
class receive_only {
  message_queue<Val>& mq;
  boost::shared_ptr<bool> alive;

  public:
  receive_only(message_queue<Val>& mq): mq(mq), alive(mq.get_alive()) {}
  boost::shared_ptr<bool> get_alive() const {return alive;}
  template <typename F>
  void receive(F continuation) {mq.receive(AMPLUSPLUS_MOVE(continuation));}
  template <typename F>
  void receive_all(F continuation) {mq.receive_all(AMPLUSPLUS_MOVE(continuation));}
  void* debug_get_queue() const {return (void*)&mq;}
};

template <typename Val, typename F>
void receive_all(receive_only<Val> from, F f) {
  AMPLUSPLUS_MOVE(from).receive_all(AMPLUSPLUS_MOVE(f));
}

template <typename Val>
void copy_messages(receive_only<Val> from, message_queue<Val>& to) {
  receive_all(from, boost::bind(&message_queue<Val>::send, boost::ref(to), _1));
}

template <typename Val, typename F>
void transform_messages(receive_only<Val> from, F f, message_queue<typename boost::result_of<const F(Val)>::type>& to) {
  receive_all(from, boost::bind(&message_queue<Val>::send, boost::ref(to), boost::bind(AMPLUSPLUS_MOVE(f), _1)));
}

}

#endif // AMPLUSPLUS_MESSAGE_QUEUE_HPP
