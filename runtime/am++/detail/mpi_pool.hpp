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

#ifndef AMPLUSPLUS_DETAIL_MPI_POOL_HPP
#define AMPLUSPLUS_DETAIL_MPI_POOL_HPP

#include <boost/pool/pool.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread.hpp>
#include <mpi.h>
#include <memory>
#include <am++/detail/mpi_global_lock.hpp>
#include <am++/detail/thread_support.hpp>

namespace amplusplus {
  namespace detail {

struct mpi_user_allocator {
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  static char* malloc(std::size_t len) {
    char* result;
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Alloc_mem(len, MPI_INFO_NULL, static_cast<void*>(&result)); AMPLUSPLUS_MPI_CALL_REGION_END
    return result;
  }
  static void free(char* p) {
    // int flag = 0;
    // MPI_Finalized(&flag);
    // if (flag) {fprintf(stderr, "FREEING AFTER FINALIZED %p\n", p); abort();}
    AMPLUSPLUS_MPI_CALL_REGION_BEGIN MPI_Free_mem(static_cast<void*>(p)); AMPLUSPLUS_MPI_CALL_REGION_END
  }
};

#if 0
class mpi_pool;

struct pool_array_deleter {
  mpi_pool& pool;
  size_t sz;
  pool_array_deleter(mpi_pool& pool, size_t sz): pool(pool), sz(sz) {}
  void operator()(void* p);
};

class mpi_pool {
  boost::pool<mpi_user_allocator> pool;
  amplusplus::detail::mutex lock;

  public:
  mpi_pool(): pool(4096), lock() {}

  boost::shared_ptr<unsigned char> alloc(size_t n) {
    size_t n_in_pages = (n + 4096 - 1) / 4096;
    boost::lock_guard<amplusplus::detail::mutex> l(lock);
    unsigned char* p = static_cast<unsigned char*>(pool.ordered_malloc(n_in_pages));
    return boost::shared_ptr<unsigned char>(p, pool_array_deleter(*this, n_in_pages));
  }

  void free(void* p, size_t sz) {
    boost::lock_guard<amplusplus::detail::mutex> l(lock);
    pool.ordered_free(p, sz);
  }
};

inline void pool_array_deleter::operator()(void* p) {
  pool.free(p, sz);
}
#endif

struct array_deleter {
  void operator()(void* p) {delete[] (unsigned char*)p;}
};

class mpi_pool {
  public:
  boost::shared_ptr<char> alloc(size_t n) {return boost::shared_ptr<char>(new char[n], boost::checked_array_deleter<char>());}
};

  }
}

#endif // AMPLUSPLUS_DETAIL_MPI_POOL_HPP
