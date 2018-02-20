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

#ifdef __bgp__
// #define BLUE_GENE_P_EXTRAS
#endif

#include <boost/smart_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/noncopyable.hpp>
#include <boost/format.hpp>
#include <boost/assert.hpp>
#include <am++/make_mpi_datatype.hpp>
#include <am++/mpi_sinha_kale_ramkumar_termination_detector.hpp>
#include <am++/mpi_sinha_kale_ramkumar_termination_detector_bgp.hpp>
#include <mpi.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <am++/mpi_transport.hpp>
#include <stdio.h>
#ifdef BLUE_GENE_P_EXTRAS
#warning BG/P mode enabled
#include <dcmf.h>
#include <dcmf_collectives.h>
#endif

#include <chrono>

// #define COLLECT_SIZE_STATS

boost::recursive_mutex amplusplus::detail::mpi_lock;

namespace amplusplus {

double get_time() {return MPI_Wtime();}

namespace detail {

#ifdef BLUE_GENE_P_EXTRAS
  struct dcmf_progress_thread {
    // From http://lists.anl-external.org/pipermail/dcmf/2010-April/000588.html
    // with modifications suggested by Jeff Hammond
    boost::shared_ptr<bool> alive;
    explicit dcmf_progress_thread(boost::shared_ptr<bool> alive): alive(AMPLUSPLUS_MOVE(alive)) {}
    void operator()() const {
      BOOST_ASSERT (alive); // Only checks that pointer is valid
      DCMF_CriticalSection_enter(0);
      while (*alive) {
        DCMF_Messager_advance();
        DCMF_CriticalSection_exit(0);
        usleep(1000);
        DCMF_CriticalSection_enter(0);
        // DCMF_CriticalSection_cycle(0);
      }
      DCMF_CriticalSection_exit(0);
    }
  };
#endif

  mpi_environment_obj::mpi_environment_obj(int argc, char ** argv, const bool need_threading, const unsigned int recv_depth, const unsigned int poll_tasks, const unsigned int flow_control_count): environment_base(), alive(boost::make_shared<bool>(true)),  recv_depth(recv_depth) , poll_tasks(poll_tasks), flow_control_count(flow_control_count) {
    int flag;
    MPI_Initialized(&flag);
    need_to_finalize_mpi = (flag == 0); // Not initialized
#ifdef BLUE_GENE_P_EXTRAS
    need_threading = true; // Since we use a progress thread
#endif
    if (need_threading) {
#ifdef AMPLUSPLUS_SINGLE_THREADED
      MPI_Init(&argc, &argv);
      std::cerr << "AMPLUSPLUS_SINGLE_THREADED defined but threading requested in mpi_environment" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
#ifdef AMPLUSPLUS_USE_THREAD_SERIALIZED
      int thrlevel;
      if (need_to_finalize_mpi) {
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thrlevel);
      } else {
	MPI_Query_thread(&thrlevel);
      }
      if (thrlevel < MPI_THREAD_SERIALIZED) {
	std::cerr << "Need MPI_THREAD_SERIALIZED or above" << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
#else
      int thrlevel;
      if (need_to_finalize_mpi) {
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thrlevel);
      } else {
	MPI_Query_thread(&thrlevel);
      }
      fprintf(stderr, "Starting with MPI_THREAD_MULTIPLE.\n");
      if (thrlevel < MPI_THREAD_MULTIPLE) {
	std::cerr << "Need MPI_THREAD_MULTIPLE or above" << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
#endif
#endif
    } else {
      if (need_to_finalize_mpi) {
	MPI_Init(&argc, &argv);
      }
    }
    {
      int thrlevel;
      MPI_Query_thread(&thrlevel);
      fprintf(stderr, "Thread level: %d\n", thrlevel);
    }
    register_builtin_mpi_datatypes();
#ifdef BLUE_GENE_P_EXTRAS
    DCMF_Messager_initialize();
    DCMF_Collective_initialize();
#endif
#ifdef BLUE_GENE_P_EXTRAS
    fprintf(stderr, "Starting BG/P progress thread.\n");
    boost::thread progress = boost::thread(dcmf_progress_thread(alive));
    progress.detach();
#endif
  }

  mpi_environment_obj::~mpi_environment_obj() {
    *alive = false;
    clear_mpi_datatype_map(); // Must be done before MPI_Finalize
    if (need_to_finalize_mpi) {
      MPI_Finalize();
#ifdef BLUE_GENE_P_EXTRAS
      DCMF_Messager_finalize();
#endif
    }
  }

  transport mpi_environment_obj::create_transport(environment& we) {
    transport t(boost::make_shared<mpi_transport_event_driven>(boost::ref(we), MPI_COMM_WORLD, recv_depth,  poll_tasks, flow_control_count), we);
#ifdef BLUE_GENE_P_EXTRAS
    t.set_termination_detector(make_mpi_sinha_kale_ramkumar_termination_detector_bgp(t));
#else
    t.set_termination_detector(make_mpi_sinha_kale_ramkumar_termination_detector(t));
#endif
    return t;
  }

  void mpi_environment_obj::set_poll_tasks(const unsigned int p) { poll_tasks = p; }
  void mpi_environment_obj::set_recv_depth(const unsigned int r) { recv_depth = r; }
  void mpi_environment_obj::set_flow_control_count(const unsigned int f) { flow_control_count = f; }
}

environment mpi_environment(int argc, char ** argv, const bool need_threading, const unsigned int recv_depth, const unsigned int poll_tasks, const unsigned int flow_control_count) {
  return environment(boost::shared_ptr<detail::mpi_environment_obj>(new detail::mpi_environment_obj(argc, argv, need_threading, recv_depth, poll_tasks, flow_control_count)));
}

void detail::swap(scoped_mpi_comm_dup& a, scoped_mpi_comm_dup& b) {std::swap(a.comm, b.comm);}

void mpi_transport_event_driven::initialize() {
  int rank_i, size_i;
  MPI_Comm_rank(comms[0], &rank_i);
  MPI_Comm_size(comms[0], &size_i);
  rank_ = (transport::rank_type)rank_i;
  size_ = (transport::rank_type)size_i;
  sends_pending_per_dest.reset(new detail::atomic<long>[size_]);
  for (transport::rank_type i = 0; i < size_; ++i) {
    sends_pending_per_dest[i].store(0);
  }
  receive_all(reqmgr.get_mpi_message_queue(),
              boost::bind(&mpi_transport_event_driven::handle_mpi_completion, this, _1));
}

void mpi_transport_event_driven::handle_mpi_completion(const detail::mpi_completion_message<detail::mpi_transport_request_info>& m) {
  const detail::mpi_transport_request_info& req_info = m.get_request_info_ref().user_info;
  if (req_info.req_kind == detail::mpi_transport_request_info::invalid_request) {
    BOOST_ASSERT (!"Invalid MPI completion");
  } else if (req_info.req_kind == detail::mpi_transport_request_info::receive_request) {
    req_info.msg_type->handle_recv_completion(m);
  } else if (req_info.req_kind == detail::mpi_transport_request_info::send_request) {
    req_info.msg_type->handle_send_completion(m);
  } else {
    BOOST_ASSERT (!"Invalid request kind");
    abort();
  }
}

void mpi_transport_event_driven::handle_termination_event(termination_message val, message_queue<termination_message>& tq) {
  // fprintf(stderr, "mpi_transport_event_driven::handle_termination_event(%d)\n", (int)val.is_last_thread());
  if (val.is_last_thread()) {
    // std::clog << (boost::format("%d mpi_transport_event_driven terminated last thread\n") % boost::this_thread::get_id()).str() << std::flush;
    for (size_t i = 0; i < this->message_types.size(); ++i) {
      this->message_types[i]->stop_receives(recvdepth, use_any_source);
    }
    // std::clog << (boost::format("%d mpi_transport_event_driven cleanup done\n") % boost::this_thread::get_id()).str() << std::flush;
  }
  BOOST_ASSERT (end_epoch_barrier);
  end_epoch_barrier->wait();
  BOOST_ASSERT (&tq);
  this->finish_end_epoch(); // May be overridden in subclasses
  tq.send(val);
}

scheduler::task_result mpi_message_type::start_one_receive_task(transport::rank_type source, size_t idx) {
  if (false /* trans.handler_calls_pending.load() > 100 */) {
    return scheduler::tr_idle;
  } else {
    this->start_one_receive(source, idx);
    return scheduler::tr_busy_and_finished;
  }
}

void mpi_message_type::start_one_receive(transport::rank_type source, size_t idx) {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  BOOST_ASSERT ((int)source == MPI_ANY_SOURCE || possible_sources->is_valid(source));
  boost::shared_ptr<void> recvbuf = this->alloc_recv_buffer();
  MPI_Request& request = this->receives[idx];
  //  fprintf(stderr, "Irecv(%p) from %d tag %zu\n", recvbuf.get(), int(source), size_t(message_index));
  AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Irecv(recvbuf.get(), this->max_count, this->get_datatype(), source, message_index, trans.comms[trans.current_comm], &request); AMPLUSPLUS_MPI_CALL_REGION_END
  trans.reqmgr.add(request, mpi_request_info<detail::mpi_transport_request_info>(detail::mpi_transport_request_info::make_receive_request(this, idx, AMPLUSPLUS_MOVE(recvbuf)), 0, message_index));
  //  fprintf(stderr, "Starting receive %p\n", request);
}

void mpi_message_type::start_receives(size_t recvdepth, bool use_any_source) {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  // std::clog << (boost::format("%d: start_receives(depth=%d, use_any_source=%d)\n") % boost::this_thread::get_id() % recvdepth % use_any_source).str() << std::flush;
  if (use_any_source) {
    this->receives.resize(recvdepth);
    for (size_t i = 0; i < recvdepth; ++i) {
      this->start_one_receive(MPI_ANY_SOURCE, i);
    }
  } else {
    size_t n = possible_sources->count();
    this->receives.resize(recvdepth * n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < recvdepth; ++j) {
        this->start_one_receive(possible_sources->rank_from_index(i), i * recvdepth + j);
      }
    }
  }
}

void mpi_message_type::stop_receives(size_t /*recvdepth*/, bool /*use_any_source*/) {
  boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
  for (size_t i = 0; i < this->receives.size(); ++i) {
    if (this->receives[i] != MPI_REQUEST_NULL) {
      // fprintf(stderr, "Canceling receive %p\n", (void*)(this->receives[i]));
      AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Cancel(&this->receives[i]); AMPLUSPLUS_MPI_CALL_REGION_END
      this->receives[i] = MPI_REQUEST_NULL;
    }
  }
  this->receives.clear();
  // std::clog << boost::this_thread::get_id() << ": ending stop_receives with " << trans.reqmgr.active_requests() << " reqs pending" << std::endl;
}

namespace {
  template <typename T>
  struct load_from_atomic {
    typedef T result_type;
    T operator()(const detail::atomic<T>* p) const {return p->load();}
  };
}

void mpi_message_type::send_untyped(const void* buf, size_t count, transport::rank_type dest, boost::function<void()> buf_deleter) {
  MPI_Datatype datatype = this->dt;
  int message_index = this->message_index;
  // fprintf(stderr, "send_untyped(%zu at %p to %zu comm %d tag %d)\n", count, buf, dest, int(trans.current_comm), int(message_index));
  BOOST_ASSERT (possible_dests->is_valid(dest));
  this->trans.td->message_send_starting(dest, message_index);
  this->trans.sends_pending_per_dest[dest].fetch_add(1);
  MPI_Request req;
  if (trans.use_ssend) {
    //    auto actstart = std::chrono::high_resolution_clock::now();
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Issend((void*)buf, count, datatype, dest, message_index, trans.comms[trans.current_comm], &req); AMPLUSPLUS_MPI_CALL_REGION_END
																	    //auto actend = std::chrono::high_resolution_clock::now();
    //    std::cout << "\t s[" << std::chrono::duration_cast<std::chrono::nanoseconds>(actend - actstart).count() << "]" << std::endl;
  } else {
    //    auto actstart = std::chrono::high_resolution_clock::now();
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Isend((void*)buf, count, datatype, dest, message_index, trans.comms[trans.current_comm], &req); AMPLUSPLUS_MPI_CALL_REGION_END
    //    auto actend = std::chrono::high_resolution_clock::now();
    //    std::cout << "\t [" << std::chrono::duration_cast<std::chrono::nanoseconds>(actend - actstart).count() << "]" << std::endl;
  }
  // fprintf(stderr, "Starting send %p\n", (void*)req);
  trans.reqmgr.add(req, mpi_request_info<detail::mpi_transport_request_info>(detail::mpi_transport_request_info::make_send_request(this, AMPLUSPLUS_MOVE(buf_deleter)), dest, message_index));
  if (trans.sends_pending_per_dest[dest].load() >= trans.flow_control_count) {
    //fprintf(stderr, "Running flow control starting at %ld requests\n", trans.sends_pending_per_dest[dest].load());
    trans.env.get_scheduler().run_until_for_flow_control(boost::bind(load_from_atomic<long>(), &this->trans.sends_pending_per_dest[dest]) < trans.flow_control_count);
  }
}

void mpi_message_type::handle_recv_completion(const detail::mpi_completion_message<detail::mpi_transport_request_info>& m) {
  const MPI_Status& st = m.get_status();
  int flag;
  AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Test_cancelled((MPI_Status*)&st, &flag); AMPLUSPLUS_MPI_CALL_REGION_END
  // fprintf(stderr, "Completed unknown receive, cancelled = %d\n", flag);
  if (flag) return;

  const mpi_request_info<detail::mpi_transport_request_info>& ri = m.get_request_info_ref();
  boost::shared_ptr<void> buf(ri.user_info.recvbuf);
  {
    boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
    this->receives[ri.user_info.receive_number] = MPI_REQUEST_NULL;
  }
  BOOST_ASSERT (st.MPI_TAG == this->message_index);
  int count;
  AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Get_count((MPI_Status*)&st, this->get_datatype(), &count); AMPLUSPLUS_MPI_CALL_REGION_END
  BOOST_ASSERT (possible_sources->is_valid(st.MPI_SOURCE));
  // fprintf(stderr, "handle_recv_completion restarted for source %d tag %d\n", st.MPI_SOURCE, st.MPI_TAG);
  trans.td->message_received(st.MPI_SOURCE, st.MPI_TAG);
  ++trans.handler_calls_pending;
  ++trans.handler_calls_pending_or_active;
  if (false /* trans.handler_calls_pending.load() > 100 */) {
    trans.env.get_scheduler().add_runnable(
      boost::bind(&mpi_message_type::start_one_receive_task,
                  this,
                  (trans.use_any_source ? MPI_ANY_SOURCE : st.MPI_SOURCE),
                  ri.user_info.receive_number));
  } else {
    this->start_one_receive((trans.use_any_source ? MPI_ANY_SOURCE : st.MPI_SOURCE), ri.user_info.receive_number);
  }
  // This needs to spawn a new task since we are in an inconsistent state inside the request manager at this point.
  handler(st.MPI_SOURCE, buf, count);
  // fprintf(stderr, "handle_recv_completion done for source %d tag %d\n", st.MPI_SOURCE, st.MPI_TAG);
}

void mpi_message_type::handle_send_completion(const detail::mpi_completion_message<detail::mpi_transport_request_info>& m) {
  const mpi_request_info<detail::mpi_transport_request_info>& ri = m.get_request_info_ref();
  ri.user_info.send_deleter();
  long cnt = this->trans.sends_pending_per_dest[ri.dest].fetch_sub(1);
  // fprintf(stderr, "Completed unknown send, pending = %ld\n", cnt);
  (void)cnt; // In case printf is commented out
  this->trans.td->message_sent(ri.dest, ri.tag);
}

bool mpi_transport_event_driven::begin_epoch() {
  // std::clog << (boost::format("%d: mpi_transport_event_driven::begin_epoch %d\n") % boost::this_thread::get_id() % this->message_types.size()).str() << std::flush;
  if (!term_queue.get()) term_queue.reset(new message_queue<termination_message>(env.get_scheduler()));
  {
    int flag;
    MPI_Status st;
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comms[(current_comm + 2) % 3], &flag, &st); AMPLUSPLUS_MPI_CALL_REGION_END
    if (flag) {
      std::cerr << (boost::format("Found mis-timed message (for previous epoch) from %d tag %d\n") % st.MPI_SOURCE % st.MPI_TAG).str() << std::flush;
      AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Abort(comms[(current_comm + 2) % 3], 9); AMPLUSPLUS_MPI_CALL_REGION_END
    }
  }
  {
    bool first = td->begin_epoch();
    if (first) {
      boost::lock_guard<amplusplus::detail::recursive_mutex> l(lock);
      current_comm = (current_comm + 1) % 3;
      size_t recvdepth = this->get_recvdepth();
      bool use_any_source = this->get_use_any_source();
      for (size_t i = 0; i < this->message_types.size(); ++i) {
        this->message_types[i]->set_message_index((int)i);
        this->message_types[i]->start_receives(recvdepth, use_any_source);
      }
    }
    // fprintf(stderr, "mpi_transport_event_driven waiting for termination on %p\n", td->get_termination_queue().debug_get_queue());
    BOOST_ASSERT (begin_epoch_barrier);
    begin_epoch_barrier->wait();
    td->get_termination_queue().receive( // For this thread only
      delay(boost::bind(&mpi_transport_event_driven::handle_termination_event, this, _1, boost::ref(*term_queue)),
            env.get_scheduler()));
    return first;
  }
}

void mpi_transport_event_driven::setup_end_epoch() {
  // fprintf(stderr, "mpi_transport_event_driven setup_end_epoch\n");
  this->flush();
  td->setup_end_epoch();
}

void mpi_transport_event_driven::setup_end_epoch_with_value(uintmax_t val) {
  // fprintf(stderr, "mpi_transport_event_driven setup_end_epoch_with_value %ju\n", val);
  /* if (td->really_ending_epoch()) */ this->flush();
  td->setup_end_epoch_with_value(val);
}

void mpi_transport_event_driven::finish_end_epoch() {}

message_type_base* mpi_transport_event_driven::create_message_type(const std::type_info& ti, size_t, transport& trans) {
  return new mpi_message_type(trans, get_mpi_datatype(ti));
}

}
