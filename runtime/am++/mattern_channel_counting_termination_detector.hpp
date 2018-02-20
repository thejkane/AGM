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

#ifndef AMPLUSPLUS_MATTERN_CHANNEL_COUNTING_TERMINATION_DETECTOR_HPP
#define AMPLUSPLUS_MATTERN_CHANNEL_COUNTING_TERMINATION_DETECTOR_HPP

#include <am++/tdetect.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread.hpp>
#include <mpi.h>

// THIS CODE IS CURRENTLY BROKEN (DEADLOCKS IN SOME CASES).

namespace amplusplus {

// Unbounded depth
// Channel counting algorithm from page 170 of
// http://www.vs.inf.ethz.ch/publ/papers/mattern-dc-1987.pdf, modified for
// non-diffusing computations

class mattern_channel_counting_termination_detector {
  boost::scoped_ptr<termination_detector> d;
  unsigned int nthreads_in_use;
  unsigned long local_value;
  mutable boost::recursive_mutex my_lock;

  public:
  template <typename Transport>
  mattern_channel_counting_termination_detector(Transport& tr) {}

  template <typename Transport>
  void initialize(Transport& tr) {
    BOOST_STATIC_ASSERT((boost::is_same<Transport, mpi_engine<mattern_channel_counting_termination_detector> >::value));
    boost::lock_guard<boost::recursive_mutex> lock(my_lock);
    nthreads_in_use = 0;
    local_value = 0;
    d.reset(new termination_detector(get_mpi_communicator(tr)));
  }

  void reset(size_t nlevels, unsigned int nthreads) {
    boost::lock_guard<boost::recursive_mutex> lock(my_lock);
    ++nthreads_in_use;
    d->reset_state();
  }

  template <typename Transport>
  void finish(Transport& tr) {
    BOOST_STATIC_ASSERT((boost::is_same<Transport, mpi_engine<mattern_channel_counting_termination_detector> >::value));
    boost::lock_guard<boost::recursive_mutex> lock(my_lock);
    if (--nthreads_in_use == 0) { // Last thread
      d->start_waiting_for_termination();
    }
    while (!d->terminated()) {
      for (size_t i = 0; i < 1000; ++i) {
        progress(tr);
        d->progress();
      }
      flush_all(tr);
    }
  }
  template <typename Transport>
  void progress_td(Transport&) {
    boost::lock_guard<boost::recursive_mutex> lock(my_lock);
    d->progress();
  }
  template <typename Transport>
  unsigned long finish(Transport& tr, unsigned long value) {
    {
      boost::lock_guard<boost::recursive_mutex> lock(my_lock);
      local_value += value;
    }
    this->finish(tr);
    MPI_Allreduce(&local_value, MPI_IN_PLACE, 1, MPI_UNSIGNED_LONG, MPI_SUM, get_mpi_communicator(tr));
    return local_value;
  }
  template <typename MsgIndexType, typename RankType>
  void before_send(MsgIndexType idx, RankType dest) {
    boost::lock_guard<boost::recursive_mutex> lock(my_lock);
    d->sending_basic_message(dest);
  }
  template <typename MsgIndexType, typename RankType>
  void message_received(MsgIndexType, RankType) {}
  template <typename MsgIndexType, typename RankType>
  void after_handler(MsgIndexType idx, RankType src) {
    boost::lock_guard<boost::recursive_mutex> lock(my_lock);
    d->finished_handling_basic_message(src);
  }
};

}

#endif // AMPLUSPLUS_MATTERN_CHANNEL_COUNTING_TERMINATION_DETECTOR_HPP
