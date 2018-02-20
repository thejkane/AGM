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
#include <mpi.h>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/ref.hpp>
#include <boost/assert.hpp>
#include <nbc.h>
#include <stdio.h>
#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <cassert>
#include <am++/mpi_transport.hpp>
#include <am++/detail/mpi_global_lock.hpp>
#include <am++/mpi_sinha_kale_ramkumar_termination_detector.hpp>
#include <iostream>

namespace amplusplus {

// Unbounded depth
// From http://charm.cs.uiuc.edu/papers/QuiescenceINTL94.pdf

// Note: initialize() is not thread-safe; everything else is
namespace {
class mpi_sinha_kale_ramkumar_termination_detector: public termination_detector_base {
  bool terminated;
  bool in_td;
  unsigned long last_total;
  enum {np_idx = 0, nc_idx = 1, user_value_idx = 2}; // Indices into *_counts
  amplusplus::detail::atomic<unsigned long> local_counts[3]; // np, nc, user value
  unsigned long global_counts[3]; // same fields
  amplusplus::detail::atomic<unsigned long> handler_starts; // For debugging
  MPI_Comm comm;
  bool start_iallreduce, iallreduce_active;
  int allreduce_start_count;
  NBC_Request reduce_req;
  int phase; // 1 or 2
  unsigned long prev_nc;
  unsigned long local_counts_to_send[3]; // copy to prevent modification of send buffer
  message_queue<termination_message> term_queue;
  mutable detail::mutex lock;
  scheduler& sched;
  transport trans;

  public:
  explicit mpi_sinha_kale_ramkumar_termination_detector(transport& _trans): term_queue(_trans.get_scheduler()), sched(_trans.get_scheduler()), trans(_trans) {
    this->initialize(trans.downcast_to_impl<mpi_transport_event_driven>()->get_mpi_communicator());
  }

  ~mpi_sinha_kale_ramkumar_termination_detector() {
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Comm_free(&comm); AMPLUSPLUS_MPI_CALL_REGION_END
  }

  receive_only<termination_message> get_termination_queue() {return term_queue;}

  private:
  void initialize(MPI_Comm comm_) {
    local_counts[user_value_idx].store(0);
    handler_starts.store(0);
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Comm_dup(comm_, &comm); AMPLUSPLUS_MPI_CALL_REGION_END
    // Initialize NBC on this comm so the allreduce won't block later
    NBC_Request req;
    {int errcode = MPI_SUCCESS; AMPLUSPLUS_NBC_CALL_REGION errcode = NBC_Ibarrier(comm, &req); if (errcode != MPI_SUCCESS) MPI_Comm_call_errhandler(this->comm, errcode);}
    {int errcode = MPI_SUCCESS; AMPLUSPLUS_NBC_CALL_REGION errcode = NBC_Wait(&req, MPI_STATUS_IGNORE); if (errcode != MPI_SUCCESS) MPI_Comm_call_errhandler(this->comm, errcode);}
    allreduce_start_count = 0;
    terminated = true;
    in_td = false;
    local_counts[np_idx].store(0);
    local_counts[nc_idx].store(0);
  }

  bool begin_epoch();
  void setup_end_epoch();
  void setup_end_epoch_with_value(uintmax_t value);
  void increase_activity_count(unsigned long);
  void decrease_activity_count(unsigned long);

  bool really_ending_epoch() const {return in_td;}

  scheduler::task_result poll_for_events(scheduler&);

  void message_being_built(size_t /*dest*/, size_t /*idx*/) {
    BOOST_ASSERT (!terminated);
    // std::cerr << "send_td -> " << dest << " tag " << idx << std::endl;
    local_counts[nc_idx].fetch_add(1);
  }
  void message_send_starting(size_t /*dest*/, size_t /*idx*/) {
  }
  void message_sent(size_t /*dest*/, size_t /*idx*/) {
  }
  void message_received(size_t /*src*/, size_t /*idx*/) {
    BOOST_ASSERT (!terminated);
    handler_starts.fetch_add(1);
  }
  void message_handled(size_t /*src*/, size_t /*idx*/) {
    BOOST_ASSERT (!terminated);
    // std::cerr << "recvd_td from " << src << " tag " << idx << std::endl;
    local_counts[np_idx].fetch_add(1);
  }
};
}
bool mpi_sinha_kale_ramkumar_termination_detector::begin_epoch() {
  // Outer code must do a thread barrier after the end of this, and only run
  // this code in one thread
  std::clog << (boost::format("%d starting mpi_sinha_kale_ramkumar_termination_detector::begin_epoch()\n") % boost::this_thread::get_id()).str() << std::flush;
  BOOST_ASSERT (terminated);
  BOOST_ASSERT (local_counts[np_idx].load() == 0);
  BOOST_ASSERT (local_counts[nc_idx].load() == 0);
  terminated = false;
  in_td = false;
  local_counts[np_idx].store(0);
  local_counts[nc_idx].store(0);
  local_counts[user_value_idx].store(0);
  handler_starts.store(0);
  std::clog << boost::this_thread::get_id() << " finishing mpi_sinha_kale_ramkumar_termination_detector::begin_epoch()" << std::endl;
  return true;
}

void mpi_sinha_kale_ramkumar_termination_detector::setup_end_epoch() {
  this->setup_end_epoch_with_value(0);
}

void mpi_sinha_kale_ramkumar_termination_detector::setup_end_epoch_with_value(uintmax_t value) {
  fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector::setup_end_epoch_with_value() top\n");
  boost::lock_guard<detail::mutex> l(this->lock);
  BOOST_ASSERT (!terminated);
  local_counts[user_value_idx].fetch_add(value);
  global_counts[np_idx] = global_counts[nc_idx] = 0;
  start_iallreduce = true;
  iallreduce_active = false;
  phase = 1;
  prev_nc = 0;
  last_total = (unsigned long)(-1);
  in_td = true;
  sched.add_idle_task(boost::bind(&mpi_sinha_kale_ramkumar_termination_detector::poll_for_events, this, _1));
  fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector::setup_end_epoch_with_value() bottom\n");
}

void mpi_sinha_kale_ramkumar_termination_detector::increase_activity_count(unsigned long v) {
  fprintf(stderr, "%p Increasing %lu\n", this, v);
  local_counts[nc_idx] += v;
}

void mpi_sinha_kale_ramkumar_termination_detector::decrease_activity_count(unsigned long v) {
  fprintf(stderr, "%p Decreasing %lu\n", this, v);
  local_counts[np_idx] += v;
}

scheduler::task_result mpi_sinha_kale_ramkumar_termination_detector::poll_for_events(scheduler&) {
  //fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector::poll_for_events running in %p count=%d\n", (void*)pthread_self(), allreduce_start_count);
  
  if (!trans.idle())
    return scheduler::tr_idle;
  
  {
    boost::lock_guard<detail::mutex> l(this->lock);
    //    fprintf(stdout, "TD progress terminated=%d iallreduce_active=%d start_iallreduce=%d\n", (int)this->terminated, (int)this->iallreduce_active, (int)this->start_iallreduce);
    if (this->terminated || !this->in_td) {
      //      fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector rp 1\n");
      return scheduler::tr_idle;
    }
    //fprintf(stderr, "Entering: entry order = %zu, last_thread = %d\n", entry_order, (int)last_thread);
    //std::clog << "Top of td loop for phase " << this->phase << std::endl;
    // std::clog << "Local side is idle (np = " << this->local_counts[this->np_idx] << ", nc = " << this->local_counts[this->nc_idx] << ", started = " << this->handler_starts << "), start_iallreduce = " << this->start_iallreduce << ", iallreduce_active = " << this->iallreduce_active << "\n" << std::flush;
#if 0
    if (0) {
      int flag;
      MPI_Status status;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, tr.get_mpi_communicator(), &flag, &status);
      std::clog << "flag = " << flag << '\n' << std::flush;
      if (flag) std::clog << "Got msg from " << status.MPI_SOURCE << " tag " << status.MPI_TAG << '\n' << std::flush;
    }
#endif
    if (this->start_iallreduce /* && this->local_counts[0].load() + this->local_counts[1].load() != last_total */ ) {
      this->prev_nc = this->global_counts[this->nc_idx];
      for (int i = 0; i < 3; ++i) this->local_counts_to_send[i] = this->local_counts[i].load();
      last_total = this->local_counts_to_send[0] + this->local_counts_to_send[1];
      //      fprintf(stderr, "Iallreduce %p %p %zu\n", (void*)this->local_counts_to_send, (void*)this->global_counts, (size_t)this->local_counts_to_send[0]);
      {int errcode = MPI_SUCCESS; AMPLUSPLUS_NBC_CALL_REGION errcode = NBC_Iallreduce((void*)this->local_counts_to_send, this->global_counts, 3, MPI_UNSIGNED_LONG, MPI_SUM, this->comm, &this->reduce_req); if (errcode != MPI_SUCCESS) MPI_Comm_call_errhandler(this->comm, errcode);}
      ++this->allreduce_start_count; // For debugging
      std::clog << (boost::format("%d: Started iallreduce phase %d np=%d nc=%d count=%d\n") % this % this->phase % this->local_counts_to_send[this->np_idx] % this->local_counts_to_send[this->nc_idx] % this->allreduce_start_count).str() << std::flush;
      this->iallreduce_active = true;
      this->start_iallreduce = false;
    }
    if (this->iallreduce_active) {
      int not_completed = NBC_OK;
      int errcode = MPI_SUCCESS;
      AMPLUSPLUS_NBC_CALL_REGION errcode = NBC_Test(&this->reduce_req, &not_completed, MPI_STATUS_IGNORE);
      if (errcode != MPI_SUCCESS) {MPI_Comm_call_errhandler(this->comm, errcode);}
      if (!not_completed) {
        std::clog << "Allreduce done A (np = " << this->global_counts[this->np_idx] << ", nc = " << this->global_counts[this->nc_idx] << ", user_value = " << this->global_counts[this->user_value_idx] << ", phase = " << this->phase << ", prev_nc = " << this->prev_nc << ")\n" << std::flush;
        this->iallreduce_active = false;
        std::clog << (boost::format("Allreduce done B (np = %d, nc = %d, user_value = %d, phase = %d, prev_nc = %d)\n") % this->global_counts[this->np_idx] % this->global_counts[this->nc_idx] % this->global_counts[this->user_value_idx] % this->phase % this->prev_nc).str() << std::flush;
        BOOST_ASSERT (this->prev_nc <= this->global_counts[this->nc_idx]); // Prevent send count from decreasing
        if (this->global_counts[this->np_idx] != this->global_counts[this->nc_idx]) {
          this->phase = 1;
          this->start_iallreduce = true;
          //          fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector rp going to phase 1\n");
          return scheduler::tr_busy;
        } else if (this->phase == 1 || this->global_counts[this->nc_idx] != this->prev_nc) {
          this->phase = 2;
          this->start_iallreduce = true;
          //          fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector rp going to phase 2\n");
          return scheduler::tr_busy;
        } else { // this->phase == 2 && this->global_counts[this->nc_idx] == this->prev_nc
          std::clog << "######Terminated\n" << std::flush;
          this->terminated = true;
          this->local_counts[this->np_idx].store(0);
          this->local_counts[this->nc_idx].store(0);
          term_queue.send(termination_message(global_counts[user_value_idx]));
          // fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector rp terminating\n");
          return scheduler::tr_remove_from_queue;
        }
        //std::clog << "Not terminated yet" << std::endl;
      } else {
        //fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector rp req incomplete\n");
        return scheduler::tr_idle;
      }
    }
    //fprintf(stderr, "mpi_sinha_kale_ramkumar_termination_detector rp retrying\n");
    return scheduler::tr_idle; // Not done yet
  }
  abort(); // Should not get here
}

termination_detector make_mpi_sinha_kale_ramkumar_termination_detector(transport& trans) {
  return boost::make_shared<mpi_sinha_kale_ramkumar_termination_detector>(boost::ref(trans));
}

}
