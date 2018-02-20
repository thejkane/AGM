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

#ifndef AMPLUSPLUS_MPI_TRANSPORT_HPP
#define AMPLUSPLUS_MPI_TRANSPORT_HPP

#include <boost/smart_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/noncopyable.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/assert.hpp>
#include <functional>
#include <memory>
#include <am++/traits.hpp>
#include <am++/make_mpi_datatype.hpp>
#include <am++/message_queue.hpp>
#include <am++/termination_detector.hpp>
#include <mpi.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <am++/detail/mpi_request_manager.hpp>
#include <am++/detail/term_detect_level_manager.hpp>
// #include <am++/detail/id_assigner.hpp>
#include <am++/detail/mpi_global_lock.hpp>
// #include <am++/detail/buffer_cache.hpp>
#include <am++/detail/thread_support.hpp>
#include <am++/transport.hpp>
#include <am++/detail/mpi_pool.hpp>
#include <am++/detail/type_info_map.hpp>

// #define COLLECT_SIZE_STATS

namespace amplusplus {

environment mpi_environment(int argc, char** argv, const bool need_threading = false, const unsigned int recvDepth = 1, const unsigned int poll_tasks = 1, const unsigned int flow_control_count = 10);

namespace detail {
  // Clone of version in Boost.MPI, with AM++ thread management
  class mpi_environment_obj : public environment_base {
  public:
    mpi_environment_obj(int argc, char ** argv, const bool need_threading, const unsigned int recv_depth, const unsigned int poll_tasks, const unsigned int flow_control_count);
    ~mpi_environment_obj();
    transport create_transport(environment& we);
    void set_poll_tasks(const unsigned int p);
    void set_recv_depth(const unsigned int r);
    void set_flow_control_count(const unsigned int f);

  private:
    bool need_to_finalize_mpi;
    boost::shared_ptr<bool> alive;
    unsigned int recv_depth;
    unsigned int poll_tasks;
    unsigned int flow_control_count;
  };
}

namespace detail {
  struct scoped_mpi_comm_dup;

  void swap(scoped_mpi_comm_dup& a, scoped_mpi_comm_dup& b);

  struct scoped_mpi_comm_dup {
    MPI_Comm comm;
    scoped_mpi_comm_dup(MPI_Comm c) {AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Comm_dup(c, &comm); AMPLUSPLUS_MPI_CALL_REGION_END}
    scoped_mpi_comm_dup(): comm(MPI_COMM_NULL) {}
#ifndef BOOST_NO_RVALUE_REFERENCES
    scoped_mpi_comm_dup(scoped_mpi_comm_dup&& o): comm(MPI_COMM_NULL) {swap(*this, o);}
#endif
    ~scoped_mpi_comm_dup() {if (comm != MPI_COMM_NULL) {AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Comm_free(&comm);AMPLUSPLUS_MPI_CALL_REGION_END}}
    operator MPI_Comm() const {return comm;}
    scoped_mpi_comm_dup& operator=(MPI_Comm c) {
      if (comm != MPI_COMM_NULL) {AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Comm_free(&comm);AMPLUSPLUS_MPI_CALL_REGION_END}
      if (c != MPI_COMM_NULL) {AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Comm_dup(c, &comm);AMPLUSPLUS_MPI_CALL_REGION_END}
      return *this;
    }

    private:
    scoped_mpi_comm_dup(const scoped_mpi_comm_dup&);
  };
}

class mpi_transport_event_driven;
class mpi_message_type;

namespace detail {
  struct mpi_transport_request_info {
    enum request_kind {invalid_request, send_request, receive_request} req_kind;
    mpi_message_type* msg_type;
    boost::function<void()> send_deleter;
    size_t receive_number;
    boost::shared_ptr<void> recvbuf;

    static mpi_transport_request_info make_send_request(mpi_message_type* msg_type_, boost::function<void()> del_) {
      mpi_transport_request_info r;
      r.req_kind = send_request;
      r.msg_type = msg_type_;
      r.send_deleter = AMPLUSPLUS_MOVE(del_);
      return r;
    }

    static mpi_transport_request_info make_receive_request(mpi_message_type* msg_type_, size_t receive_number_, const boost::shared_ptr<void>& recvbuf_) {
      mpi_transport_request_info r;
      r.req_kind = receive_request;
      r.msg_type = msg_type_;
      r.receive_number = receive_number_;
      r.recvbuf = recvbuf_;
      return r;
    }

    void swap(mpi_transport_request_info& o) {
      std::swap(req_kind, o.req_kind);
      std::swap(msg_type, o.msg_type);
      send_deleter.swap(o.send_deleter);
      std::swap(receive_number, o.receive_number);
      recvbuf.swap(o.recvbuf);
    }

    friend std::ostream& operator<<(std::ostream& o, const mpi_transport_request_info& r) {
      o << "{" << (r.req_kind == mpi_transport_request_info::send_request ? "send" : r.req_kind == mpi_transport_request_info::receive_request ? "recv" : "unknown") << "}";
      return o;
    }

    // private: mpi_transport_request_info() {}
  };
}

class mpi_transport_event_driven: public transport_base, boost::noncopyable {
  public:
  explicit mpi_transport_event_driven(environment& env, MPI_Comm comm = MPI_COMM_WORLD, int recvDepth =1 , int poll_tasks = 1, int flow_control_count = 10)
    : env(env), reqmgr(env.get_scheduler(), poll_tasks), current_comm(0),
      recvdepth(recvDepth), nthreads(1), use_any_source(false), use_ssend(false),
      begin_epoch_barrier(new boost::barrier(1)),
      end_epoch_barrier(new boost::barrier(1)),
      term_queue(), flow_control_count(flow_control_count)
  {
    BOOST_ASSERT(flow_control_count > 0);
    comms[0] = comm;
    comms[1] = comm;
    comms[2] = comm;
    this->initialize();
  }

  ~mpi_transport_event_driven() {
#ifndef NDEBUG
    for (size_t i = 0; i < size(); ++i) {
      BOOST_ASSERT (sends_pending_per_dest[i].load() == 0);
    }
#endif
  }

  private:
  void initialize();
  void handle_mpi_completion(const detail::mpi_completion_message<detail::mpi_transport_request_info>& m);
  void handle_termination_event(termination_message val, message_queue<termination_message>& tq);

  private:
  friend class mpi_message_type;

  public:
  bool begin_epoch();
  void setup_end_epoch();
  void setup_end_epoch_with_value(uintmax_t val);
  void finish_end_epoch();

  boost::shared_ptr<void> alloc_memory(size_t sz) const {
    return boost::static_pointer_cast<void>(pool.alloc(sz));
  }

  // This communicator should be MPI_Comm_dup'ed before use
  MPI_Comm get_mpi_communicator() const {
    return comms[current_comm];
  }

  bool is_valid_rank(rank_type r) const {return r < size();}

  void set_use_ssend(bool x) {use_ssend = x;}
  bool get_use_ssend() const {return use_ssend;}
  void set_recvdepth(size_t recvdepth_) {BOOST_ASSERT (recvdepth_ >= 1); recvdepth = recvdepth_;}
  size_t get_recvdepth() const {return recvdepth;}
  void set_use_any_source(bool use_any_source_) {use_any_source = use_any_source_;}
  bool get_use_any_source() const {return use_any_source;}

  void set_nthreads(size_t nt) {
    nthreads = nt;
    td->set_nthreads(nt);
    begin_epoch_barrier.reset(new boost::barrier(nt));
    end_epoch_barrier.reset(new boost::barrier(nt));
  }
  size_t get_nthreads() const {return nthreads;}

  message_type_base* create_message_type(const std::type_info& ti, size_t size, transport&);

  const transport_base::rank_type& rank() const {return rank_;}
  const transport_base::rank_type& size() const {return size_;}
  
  // bool any_sends_pending() const {return sends_pending.load();}

  void set_termination_detector(const termination_detector& td_) {
    if (boost::shared_ptr<detail::td_thread_wrapper> w = boost::dynamic_pointer_cast<detail::td_thread_wrapper>(td_)) {
      td = w;
    } else {
      td = boost::make_shared<detail::td_thread_wrapper>(td_, boost::ref(env.get_scheduler()));
    }
  }
  termination_detector get_termination_detector() const {return td;}

  receive_only<termination_message> get_termination_queue() { // Per thread
    if (!term_queue.get()) term_queue.reset(new message_queue<termination_message>(env.get_scheduler()));
    return *term_queue;
  }

  void increase_activity_count(unsigned long long v) {td->increase_activity_count(v);}
  void decrease_activity_count(unsigned long long v) {td->decrease_activity_count(v);}
  bool do_explicit_polling(){//std::cout<<"Calling poll function from mpi_transport_event_driven\n";
    return reqmgr.do_poll_explicitly();
  }

  private:
  void add_message_type(mpi_message_type* mt) {
    BOOST_ASSERT (!td->in_epoch());
    message_types.push_back(mt);
  }
  void remove_message_type(mpi_message_type* mt) {
    BOOST_ASSERT (!td->in_epoch());
    std::vector<mpi_message_type*>::iterator i = std::find(message_types.begin(), message_types.end(), mt);
    BOOST_ASSERT (i != message_types.end());
    message_types.erase(i);
  }
  void replace_message_type(mpi_message_type* old_mt, mpi_message_type* new_mt) {
    BOOST_ASSERT (!td->in_epoch());
    std::vector<mpi_message_type*>::iterator i = std::find(message_types.begin(), message_types.end(), old_mt);
    BOOST_ASSERT (i != message_types.end());
    *i = new_mt;
  }

  private:
  environment& env;
  transport::rank_type rank_, size_;
  detail::mpi_request_manager<detail::mpi_transport_request_info> reqmgr;
  int current_comm;
  detail::scoped_mpi_comm_dup comms[3];
  size_t recvdepth;
  size_t nthreads;
  bool use_any_source, use_ssend;
  boost::scoped_array<detail::atomic<long> > sends_pending_per_dest;
  std::vector<mpi_message_type*> message_types;
  amplusplus::detail::recursive_mutex lock;
  mutable detail::mpi_pool pool;
  boost::shared_ptr<detail::td_thread_wrapper> td;
  boost::scoped_ptr<boost::barrier> begin_epoch_barrier;
  boost::scoped_ptr<boost::barrier> end_epoch_barrier;
  boost::thread_specific_ptr<message_queue<termination_message> > term_queue;
  const int flow_control_count;
};

class mpi_message_type: public message_type_base, public boost::noncopyable {
  public:
  explicit mpi_message_type(transport trans, MPI_Datatype dt)
    : message_type_base(trans), valid(true), trans(*trans.downcast_to_impl<mpi_transport_event_driven>()), trans_wrapped(trans), dt(dt), handler(), max_count(0), possible_dests(), possible_sources()
  {
    this->trans.add_message_type(this);
    MPI_Aint lb, sz;
    MPI_Type_get_extent(dt, &lb, &sz);
    dt_size = (size_t)sz;
  }

  virtual ~mpi_message_type() {
    if (this->valid) {
      BOOST_ASSERT (this->receives.empty());
      trans.remove_message_type(this);
      this->valid = false;
    }
  }
  boost::shared_ptr<void> alloc_recv_buffer() const {
    return trans.alloc_memory(this->max_count * this->dt_size);
  }
  MPI_Datatype get_datatype() const {return dt;}
  // mpi_transport_event_driven& transport() const {return trans;}
  transport& get_transport() const {return trans_wrapped;}

  virtual void message_being_built(transport::rank_type dest) {trans_wrapped.message_being_built(dest, message_index);}
  virtual void handler_done(transport::rank_type src) {trans.td->message_handled(int(src), this->message_index);}
  virtual bool flush(transport::rank_type /*dest*/) {return false;}
  virtual scheduler::task_result flush_all() {return scheduler::tr_idle;}

  void start_one_receive(transport::rank_type source, size_t idx);
  scheduler::task_result start_one_receive_task(transport::rank_type source, size_t idx);
  void start_receives(size_t recvdepth, bool use_any_source);
  void stop_receives(size_t recvdepth, bool use_any_source);
  void send_untyped(const void* buf, size_t count, transport::rank_type dest, boost::function<void()> buf_deleter);

  void handle_recv_completion(const detail::mpi_completion_message<detail::mpi_transport_request_info>& m);
  void handle_send_completion(const detail::mpi_completion_message<detail::mpi_transport_request_info>& m);

  void set_handler_internal(message_type_base::handler_type h) {handler = AMPLUSPLUS_MOVE(h);}

  void set_message_index(int idx) {message_index = idx;}

  void set_max_count(size_t m) {max_count = m;}
  size_t get_max_count() const {return max_count;}
  void set_possible_sources(valid_rank_set p) {possible_sources = p;}
  valid_rank_set get_possible_sources() const {return possible_sources;}
  void set_possible_dests(valid_rank_set p) {possible_dests = p;}
  valid_rank_set get_possible_dests() const {return possible_dests;}

  protected:
  bool valid;
  mpi_transport_event_driven& trans;
  mutable transport trans_wrapped;
  MPI_Datatype dt;
  size_t dt_size;
  int message_index;
  std::vector<MPI_Request> receives;
  message_type_base::handler_type handler;
  size_t max_count;
  valid_rank_set possible_dests;
  valid_rank_set possible_sources;
  // Only within-epoch calls are locked (others are not thread-safe)
  amplusplus::detail::recursive_mutex lock;
};

}

#endif // AMPLUSPLUS_MPI_TRANSPORT_HPP
