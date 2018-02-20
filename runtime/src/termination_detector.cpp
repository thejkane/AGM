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

#include <config.h>

#include <am++/termination_detector.hpp>
#include <am++/traits.hpp>
#include <am++/message_queue.hpp>
#include <am++/detail/thread_support.hpp>
#include <boost/format.hpp>
#include <boost/assert.hpp>

namespace amplusplus {
namespace detail {

td_thread_wrapper::td_thread_wrapper(termination_detector td, scheduler& sched)
  : td(td),
    begin_epoch_barrier(new detail::barrier(1)),
    nthreads_total(1),
    msg_queues(1, (message_queue<termination_message>*)NULL),
    sched(sched)
{
  local_finish_value = 0;
  thread_id_counter.store(0);
  nthreads_active.store(0);
  nthreads_in_epoch.store(0);
  currently_in_epoch.store(0);
}

bool td_thread_wrapper::begin_epoch() {
  const int nthr = nthreads_total;

  if (!thread_data_ptr.get()) thread_data_ptr.reset(new per_thread_data(sched));
  per_thread_data& my_thread_data = *thread_data_ptr;

  my_thread_data.id = (thread_id_counter.fetch_add(1) % nthr);
  BOOST_ASSERT (my_thread_data.activity_count == 0);
  my_thread_data.activity_count = 0;
  BOOST_ASSERT (my_thread_data.id >= 0 && my_thread_data.id < nthr && my_thread_data.id < (int)msg_queues.size());
  msg_queues[my_thread_data.id] = &my_thread_data.term_queue;

  if (my_thread_data.id == 0) {
    boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
    local_finish_value = 0;
    nthreads_active.store(0);
    nthreads_in_epoch.store(nthreads_total);
    currently_in_epoch.store(1);
    td->begin_epoch();
  }

  begin_epoch_barrier->wait();
  return my_thread_data.id == 0;
}

void td_thread_wrapper::setup_end_epoch() {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  if (nthreads_in_epoch.fetch_add(-1) == 1) {
    td->setup_end_epoch();
    // fprintf(stderr, "thread_wrapper setup_end_epoch\n");
    td->get_termination_queue().receive(boost::bind(&td_thread_wrapper::handle_termination_message, this, _1));
  }
}

void td_thread_wrapper::setup_end_epoch_with_value(uintmax_t val) {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  local_finish_value += val;
  if (nthreads_in_epoch.fetch_add(-1) == 1) {
    td->setup_end_epoch_with_value(local_finish_value);
    // fprintf(stderr, "thread_wrapper setup_end_epoch_with_value\n");
    td->get_termination_queue().receive(boost::bind(&td_thread_wrapper::handle_termination_message, this, _1));
  }
}

void td_thread_wrapper::handle_termination_message(termination_message msg) {
  // fprintf(stderr, "thread_wrapper handle_termination_message\n");
  BOOST_ASSERT (msg_queues.size() == (size_t)nthreads_total);
  BOOST_ASSERT (nthreads_in_epoch.load() == 0);
  BOOST_ASSERT (msg.is_last_thread());
  currently_in_epoch.store(0);
  for (int i = 0; i < nthreads_total; ++i) {
    BOOST_ASSERT (msg_queues[i]);
    // fprintf(stderr, "Sending termination to %p\n", msg_queues[i]);
    msg_queues[i]->send(termination_message(msg.get_combined_value(), (i == nthreads_total - 1)));
  }
}

bool td_thread_wrapper::in_epoch() const {
  return currently_in_epoch.load() != 0;
}

bool td_thread_wrapper::really_ending_epoch() const {
  std::cout << "inside really_ending_epoch ..." << std::endl;
  return nthreads_active.load() == 0 &&
         nthreads_in_epoch.load() == 0 &&
         currently_in_epoch.load() == 1;
}

void td_thread_wrapper::set_nthreads(size_t n) {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  nthreads_total = (int)n;
  thread_id_counter.store(0);
  msg_queues.resize(n);
  begin_epoch_barrier.reset(new detail::barrier(n));
}

size_t td_thread_wrapper::get_nthreads() const {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  return nthreads_total;
}

}
}
