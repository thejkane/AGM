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

#ifndef AMPLUSPLUS_TERMINATION_DETECTOR_HPP
#define AMPLUSPLUS_TERMINATION_DETECTOR_HPP

#include <am++/traits.hpp>
#include <am++/detail/thread_support.hpp>
#include <am++/message_queue.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/format.hpp>
#include <boost/assert.hpp>
#include <iostream>

namespace amplusplus {

struct termination_message {
  public:
  // Default for second parameter is because non-threaded TD implementations always put true there
  termination_message(uintmax_t val, bool last_terminating_thread = true): val(val), last_terminating_thread(last_terminating_thread) {}
  uintmax_t get_combined_value() const {return val;}
  bool is_last_thread() const {return last_terminating_thread;}

  private:
  uintmax_t val;
  bool last_terminating_thread;
};

class termination_detector_base {
  public:
  virtual bool begin_epoch() = 0; // Returns true in first thread to get to begin_epoch
  virtual void setup_end_epoch() = 0;
  virtual void setup_end_epoch_with_value(uintmax_t val) = 0;
  virtual bool really_ending_epoch() const = 0;
  virtual void increase_activity_count(unsigned long) = 0; // Something local has happened that might send messages
  virtual void decrease_activity_count(unsigned long) = 0; // That has stopped
  virtual receive_only<termination_message> get_termination_queue() = 0;
  virtual void message_received(size_t source, size_t msg_type) = 0;
  virtual void message_handled(size_t source, size_t msg_type) = 0;
  virtual void message_being_built(size_t dest, size_t msg_type) = 0;
  virtual void message_send_starting(size_t dest, size_t msg_type) = 0;
  virtual void message_sent(size_t dest, size_t msg_type) = 0;
  virtual void set_nthreads(size_t n = 1) {
    BOOST_ASSERT (n == 1);
    (void)n;
  }
  virtual size_t get_nthreads() const {return 1;}
  virtual ~termination_detector_base() {}
};

typedef boost::shared_ptr<termination_detector_base> termination_detector;

namespace detail {
class td_thread_wrapper: public amplusplus::termination_detector_base {
  public:
  td_thread_wrapper(termination_detector td, scheduler&);

  public:
  bool begin_epoch();
  void setup_end_epoch();
  void setup_end_epoch_with_value(uintmax_t val);

  receive_only<termination_message> get_termination_queue() { // Returns different value in each thread
    if (!thread_data_ptr.get()) thread_data_ptr.reset(new per_thread_data(sched));
    return thread_data_ptr->term_queue;
  }

  bool in_epoch() const;
  bool really_ending_epoch() const;

  void increase_activity_count(unsigned long v) {
    if (v == 0) return;
    BOOST_ASSERT (thread_data_ptr.get());
    unsigned long long& ac = thread_data_ptr->activity_count;
    unsigned long long old_val = ac;
    ac += v;
    std::clog << (boost::format("%x increase_activity_count(%d) from %d to %d\n") % this % v % old_val % ac).str() << std::flush;
    if (old_val == 0) {if (nthreads_active.fetch_add(1) == 0) td->increase_activity_count(1);}
  }

  void decrease_activity_count(unsigned long v) {
    if (v == 0) return;
    BOOST_ASSERT (thread_data_ptr.get());
    unsigned long long& ac = thread_data_ptr->activity_count;
    BOOST_ASSERT (v <= ac);
    ac -= v;
    std::clog << (boost::format("%x decrease_activity_count(%d) from %d to %d\n") % this % v % (ac + v) % ac).str() << std::flush;
    if (ac == 0) {if (nthreads_active.fetch_add(-1) == 1) td->decrease_activity_count(1);}
  }

  void message_received(size_t source, size_t msg_type) {td->message_received(source, msg_type);}
  void message_handled(size_t source, size_t msg_type) {td->message_handled(source, msg_type);}
  void message_being_built(size_t dest, size_t msg_type) {td->message_being_built(dest, msg_type);}
  void message_send_starting(size_t dest, size_t msg_type) {td->message_send_starting(dest, msg_type);}
  void message_sent(size_t dest, size_t msg_type) {td->message_sent(dest, msg_type);}

  void set_nthreads(size_t n);
  size_t get_nthreads() const;

  private:
  termination_detector td;
  // Encoding of epoch_barrier_status (see promela/mpi_local_threads.pr for an old version of this):
  // 0..nthreads-1 = waiting for everyone to enter epoch,
  // nthreads..2*nthreads-1 = waiting for everyone to leave epoch,
  // 2*nthreads = waiting for remote ranks to finish
  // 2*nthreads+1..3*nthreads = waiting for all local threads to know about termination
  boost::shared_ptr<detail::barrier> begin_epoch_barrier;
  detail::atomic<int> nthreads_active;
  detail::atomic<int> nthreads_in_epoch;
  int nthreads_total;
  detail::atomic<int> currently_in_epoch;
  uintmax_t local_finish_value;
  detail::atomic<int> thread_id_counter;
  struct per_thread_data {
    int id;
    unsigned long long activity_count;
    message_queue<termination_message> term_queue;
    per_thread_data(scheduler& sched): id(-1), activity_count(0), term_queue(sched) {}
  };
  boost::thread_specific_ptr<per_thread_data> thread_data_ptr;
  std::vector<message_queue<termination_message>*> msg_queues; // Pointers into per-thread data to allow access by other threads
  scheduler& sched;
  mutable amplusplus::detail::recursive_mutex lock;

  void handle_termination_message(termination_message msg);
};

}

}

#endif // AMPLUSPLUS_TERMINATION_DETECTOR_HPP
