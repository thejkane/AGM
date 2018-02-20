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

#ifndef AMPLUSPLUS_DETAIL_MPI_GLOBAL_LOCK_HPP
#define AMPLUSPLUS_DETAIL_MPI_GLOBAL_LOCK_HPP

#include <boost/thread.hpp>
#include <boost/noncopyable.hpp>
#include <am++/detail/thread_support.hpp>

#ifdef AMPLUSPLUS_SINGLE_THREADED
#ifdef AMPLUSPLUS_USE_THREAD_SERIALIZED
#error "AMPLUSPLUS_SINGLE_THREADED and AMPLUSPLUS_USE_THREAD_SERIALIZED conflict"
#endif
#endif

namespace amplusplus {
  namespace detail {
    // Even under MPI_THREAD_MULTIPLE, we need to guard calls to libNBC
    extern boost::recursive_mutex mpi_lock;

    struct mpi_lock_guard_for_if {
      boost::lock_guard<boost::recursive_mutex> x;
      mpi_lock_guard_for_if(int): x(mpi_lock) {}
      mpi_lock_guard_for_if(const mpi_lock_guard_for_if&): x(mpi_lock) {}
      operator bool() const {return false;}
    };
  }
}

#ifdef AMPLUSPLUS_USE_THREAD_SERIALIZED
#define AMPLUSPLUS_MPI_CALL_REGION_BEGIN \
  { ::amplusplus::detail::mpi_lock_guard_for_if lg = 0;
#define AMPLUSPLUS_MPI_CALL_REGION_END }
#else // AMPLUSPLUS_USE_THREAD_SERIALIZED
#define AMPLUSPLUS_MPI_CALL_REGION_BEGIN {
#define AMPLUSPLUS_MPI_CALL_REGION_END }
#endif
#define AMPLUSPLUS_NBC_CALL_REGION					\
  if (::amplusplus::detail::mpi_lock_guard_for_if lg = 0) {} else

#endif // AMPLUSPLUS_DETAIL_MPI_GLOBAL_LOCK_HPP
