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

#ifndef AMPLUSPLUS_DETAIL_MPI_REQUEST_MANAGER_HPP
#define AMPLUSPLUS_DETAIL_MPI_REQUEST_MANAGER_HPP

#include <am++/traits.hpp>
#include <am++/message_queue.hpp>
#include <am++/detail/thread_support.hpp>
#include <am++/detail/mpi_global_lock.hpp>
#include <vector>
#include <list>
#include <utility>
#include <typeinfo>
#include <mpi.h>
#include <boost/thread.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/assert.hpp>

namespace amplusplus {

template <typename UserInfo>
struct mpi_request_info {
  bool slot_has_active_request; // Used instead of testing for MPI_REQUEST_NULL to allow incremental sending to output queue of multiple completed requests without things getting overwritten in the request_info buffer
  size_t dest;
  size_t tag;
  UserInfo user_info;
  mpi_request_info(): slot_has_active_request(false), dest(size_t(-1)), tag(size_t(-1)), user_info() {}
  mpi_request_info(UserInfo user_info, size_t dest, size_t tag): slot_has_active_request(true), dest(dest), tag(tag), user_info(AMPLUSPLUS_MOVE(user_info)) {}
  friend std::ostream& operator<<(std::ostream& o, const mpi_request_info& r) {
    o << (r.slot_has_active_request ? "active" : "inactive") << " dest " << r.dest << " tag " << r.tag << " user_info " << r.user_info;
    return o;
  }
  void swap(mpi_request_info& o) {
    std::swap(slot_has_active_request, o.slot_has_active_request);
    std::swap(dest, o.dest);
    std::swap(tag, o.tag);
    user_info.swap(o.user_info);
  }
};

  namespace detail {

template <typename UserInfo>
struct mpi_completion_message {
  public:
  explicit mpi_completion_message(const MPI_Status& status = MPI_Status(), mpi_request_info<UserInfo> info = mpi_request_info<UserInfo>()): info(AMPLUSPLUS_MOVE(info)), status(status) {}
  MPI_Status get_status() const {return status;}
  const mpi_request_info<UserInfo>& get_request_info_ref() const {return info;}
  mpi_request_info<UserInfo> get_request_info() const {return info;}

  private:
  mpi_request_info<UserInfo> info;
  MPI_Status status;
};

template <typename UserInfo>
class mpi_request_manager: boost::noncopyable {
  private:
  std::vector<MPI_Request> reqs;
  std::vector<mpi_request_info<UserInfo> > req_info;
  std::vector<int> indices;
  std::vector<MPI_Status> statuses;
  int nreqs_active;
  message_queue<mpi_completion_message<UserInfo> > mpi_msg_queue;
  atomic<int> add_pending;
  boost::shared_ptr<bool> need_to_exit;
  amplusplus::scheduler& sched;

  mutable amplusplus::detail::recursive_mutex req_queue_lock; // External so it can share an external MPI lock

  scheduler::task_result poll_for_messages(scheduler& sched, boost::shared_ptr<bool> need_to_exit);
  scheduler::task_result poll_unlocked(boost::unique_lock<recursive_mutex>&);

  public:
  mpi_request_manager(scheduler& sched, int poll_tasks): nreqs_active(0), mpi_msg_queue(sched), need_to_exit(new bool(false)), sched(sched) {
    BOOST_ASSERT(poll_tasks > 0);
    add_pending.store(0);
    for(int i = 0; i != poll_tasks; ++i)
      sched.add_runnable(boost::bind(&mpi_request_manager::poll_for_messages, this, _1, need_to_exit)); // Extra argument needed to force a copy
  }
  ~mpi_request_manager() {/* fprintf(stderr, "~mpi_request_manager() on %p\n", this); */ *need_to_exit = true;}

  void add(MPI_Request req, mpi_request_info<UserInfo> info);
  receive_only<mpi_completion_message<UserInfo> > get_mpi_message_queue() {return mpi_msg_queue;}

  bool empty() const {
    boost::lock_guard<amplusplus::detail::recursive_mutex> lock(req_queue_lock);
    return nreqs_active == 0 && mpi_msg_queue.empty() && add_pending.load() == 0;
  }

  int size() const {
    boost::lock_guard<amplusplus::detail::recursive_mutex> lock(req_queue_lock);
    return nreqs_active;
  }

  bool do_poll_explicitly()
  {
    return poll_for_messages(sched,need_to_exit) == scheduler::tr_busy;
  }
};

template <typename UserInfo>
void mpi_request_manager<UserInfo>::add(MPI_Request req, mpi_request_info<UserInfo> info) {
  if (req == MPI_REQUEST_NULL) return;
  // std::clog << (boost::format("%x add %d with %d\n") % uintptr_t(this) % uintptr_t(req) % info).str() << std::flush;
  ++add_pending;
  {
    boost::lock_guard<amplusplus::detail::recursive_mutex> lock(req_queue_lock);
    ++nreqs_active;
    --add_pending;
    for (size_t i = 0; i < reqs.size(); ++i) {
      if (!req_info[i].slot_has_active_request) {
        reqs[i] = req;
        req_info[i] = AMPLUSPLUS_MOVE(info);
        req_info[i].slot_has_active_request = true;
        // std::clog << (boost::format("%x add %d with %d found slot %d\n") % uintptr_t(this) % uintptr_t(req) % info % i).str() << std::flush;
        goto do_poll;
      }
    }
    reqs.push_back(req);
    req_info.push_back(AMPLUSPLUS_MOVE(info));
    req_info.back().slot_has_active_request = true;
    // std::clog << (boost::format("%x add %d with %d allocating slot %d\n") % uintptr_t(this) % uintptr_t(req) % info % (reqs.size() - 1)).str() << std::flush;
    indices.push_back(0);
    statuses.push_back(MPI_Status());
    // std::clog << "add " << req << " -> " << reqs.size() << '\n' << std::flush;
    do_poll: ; // this->poll_unlocked();
    // fprintf(stderr, "Add on %p updated count to %zu of %zu\n", this, size_t(nreqs_active), size_t(reqs.size()));
  }
  // sched.run_all_tail();
}

template <typename UserInfo>
scheduler::task_result mpi_request_manager<UserInfo>::poll_for_messages(scheduler&, boost::shared_ptr<bool> need_to_exit) {
  // fprintf(stderr, "mpi_request_manager<%s, %p>::poll_for_messages trying to run\n", typeid(UserInfo).name(), this);
  if (*need_to_exit) {/* fprintf(stderr, "poll_for_messages %p removing itself from queue\n", this); */ return scheduler::tr_remove_from_queue;}
  if (add_pending.load() != 0) return scheduler::tr_idle; // Don't try to lock in competition with an add
  {
    boost::unique_lock<amplusplus::detail::recursive_mutex> lock(req_queue_lock, boost::try_to_lock);
    if (!lock.owns_lock()) return scheduler::tr_idle;
    return this->poll_unlocked(lock);
  }
}

namespace detail {
  template <typename UserInfo>
  struct send_range_task {
    message_queue<mpi_completion_message<UserInfo> >& mpi_msg_queue;
    std::vector<mpi_completion_message<UserInfo> > msgs;
    send_range_task(message_queue<mpi_completion_message<UserInfo> >& mpi_msg_queue, const std::vector<mpi_completion_message<UserInfo> >& msgs): mpi_msg_queue(mpi_msg_queue), msgs(msgs) {}
    scheduler::task_result operator()(scheduler&) const {
      mpi_msg_queue.send_range(msgs.begin(), msgs.end());
      return scheduler::tr_busy_and_finished;
    }
  };
}

template <typename UserInfo>
scheduler::task_result
mpi_request_manager<UserInfo>::poll_unlocked(boost::unique_lock<amplusplus::detail::recursive_mutex>& l) {
  // fprintf(stderr, "mpi_request_manager<%s, %p>::poll_unlocked running %zu of %zu\n", typeid(UserInfo).name(), this, size_t(nreqs_active), size_t(reqs.size()));
  assert (l.owns_lock()); (void)l;
  if (nreqs_active == 0) return scheduler::tr_idle;
  int outcount;
  if (!reqs.empty()) { // Avoid Open MPI failure when output arrays are NULL
    BOOST_ASSERT (indices.size() >= reqs.size());
    BOOST_ASSERT (statuses.size() >= reqs.size());
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Testsome((int)reqs.size(), &reqs[0], &outcount, &indices[0], &statuses[0]); AMPLUSPLUS_MPI_CALL_REGION_END
    // fprintf(stderr, "MPI_Testsome returned %d\n", outcount);
  } else {
    BOOST_ASSERT (this->nreqs_active == 0);
    return scheduler::tr_idle;
  }
  if (outcount == 0) {
    return scheduler::tr_idle;
  } else if (outcount == MPI_UNDEFINED) {
    BOOST_ASSERT (this->nreqs_active == 0);
    return scheduler::tr_idle;
  }
  nreqs_active -= outcount;
  for (int i = 0; i < outcount; ++i) {
    mpi_request_info<UserInfo> ri;
    ri.swap(req_info[indices[i]]);
    // MPI message queue receiver cannot directly send messages (and thus
    // shouldn't call handlers) without spawning a task.
    mpi_msg_queue.send(mpi_completion_message<UserInfo>(statuses[i], AMPLUSPLUS_MOVE(ri)));
  }
  return scheduler::tr_busy;
}

  }
}

#endif // AMPLUSPLUS_DETAIL_MPI_REQUEST_MANAGER_HPP
