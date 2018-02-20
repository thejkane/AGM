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

#ifndef AMPLUSPLUS_DETAIL_BUFFER_CACHE_HPP
#define AMPLUSPLUS_DETAIL_BUFFER_CACHE_HPP

#include <boost/smart_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/assert.hpp>
#include <iostream>
#include <memory>
#include <boost/noncopyable.hpp>
#include <am++/detail/thread_support.hpp>

namespace amplusplus {
  namespace detail {

#if 0
class buffer_cache: boost::noncopyable {
  std::vector<boost::shared_ptr<void> > all_buffers; // Keeps ownership of all of them
  boost::mutex all_buffers_lock;
  // Entries in buffer_ptrs are 0 for "unallocated", (void*)1 for "waiting for
  // someone else to do allocation", and a pointer for "allocated"
  amplusplus::detail::atomic<amplusplus::detail::atomic<void*>*> buffer_ptrs[20];
  amplusplus::detail::atomic<size_t> last_entry_written_plus_1;
  mpi_pool& pool;

  static const size_t page_0_size = (1 << 15);
  static const int num_pages = 20;

  static size_t page_size(size_t i) {
    return page_0_size << ((i == 0) ? 0 : (i - 1));
  }

  static size_t page_start(size_t i) {
    return (i == 0) ? 0 : (page_0_size << (i - 1));
  }

  struct free_buf {
    amplusplus::detail::atomic<void*>* loc_in_page_list;
    free_buf(amplusplus::detail::atomic<void*>* loc_in_page_list): loc_in_page_list(loc_in_page_list) {}
    void operator()(void* b) {
      void* old_val = loc_in_page_list->swap(b);
      BOOST_ASSERT (old_val == 0);
    }
  };

  public:
  buffer_cache(mpi_pool& pool): last_entry_written_plus_1(0), pool(pool) {
    buffer_ptrs[0].store(new amplusplus::detail::atomic<void*>[page_0_size]);
    for (int i = 1; i < num_pages; ++i) buffer_ptrs.store(0);
  }

  ~buffer_cache() {
    // The buffers themselves are owned by all_buffers so we don't need to free
    // them explicitly
    for (int i = 0; i < num_pages; ++i) delete[] buffer_ptrs[i].load();
  }

  boost::shared_ptr<void> allocate() {
    void* p = 0;
    amplusplus::detail::atomic<void*>* p_ptr = 0;
    size_t index = 0;
    for (int i = 0; i < num_pages; ++i) {
      amplusplus::detail::atomic<void*>* page = 0;
      const size_t ps = page_size(i);
      while (true) {
        const bool changed = buffer_ptrs[i].compare_exchange_strong(page, (amplusplus::detail::atomic<void*>*)1);
        if (changed) {
          page = new amplusplus::detail::atomic<void*>[ps];
          for (size_t i = 0; i < ps; ++i) page[i].store(0);
          buffer_ptrs[i].store(page);
          break;
        }
        if (page == (amplusplus::detail::atomic<void*>*)1) continue; // Spin until allocation by another thread is done
        break;
      }
      bool last_page = (last_entry_written_plus_1.load() - page_start(i) < ps);
      const size_t j_end = (last_page ? last_entry_written_plus_1.load() - page_start(i) : ps);
      for (size_t j = 0; j < j_end; ++j) {
        p = page[j].swap(0);
        if (p != 0) {
          p_ptr = &page[j];
          goto have_buffer;
        }
      }
      if (last_page) {
        p_ptr = ...
      }
    }
    need_new_buffer: {
      boost::shared_ptr<void> buf = pool.alloc();
      {
        boost::scoped_lock<boost::mutex> l(all_buffers_lock);
        all_buffers.push_back(buf);
      }
      return buf.get();
    }
    have_buffer: {
      return boost::shared_ptr<void>(p, free_buf(p_ptr));
    }
  }
};
#endif

#if 1
class buffer_cache: boost::noncopyable {
  static const size_t size = (1 << 16);

  std::vector<boost::shared_ptr<void> > all_buffers; // Keeps ownership of all of them
  amplusplus::detail::mutex all_buffers_lock;
  amplusplus::detail::atomic<void*> buffer_ptrs[size];
  amplusplus::detail::atomic<size_t> last_entry_written_plus_1;
  transport trans;
  const size_t buffer_size;
  boost::shared_ptr<bool> cache_deleted;

  struct free_buf {
    boost::shared_ptr<bool> cache_deleted;
    amplusplus::detail::atomic<void*>* p_ptr;
    free_buf(boost::shared_ptr<bool> cache_deleted, amplusplus::detail::atomic<void*>* p_ptr)
      : cache_deleted(cache_deleted), p_ptr(p_ptr) {}
    void operator()(void* b) {
      if (*cache_deleted) return;
#ifndef NDEBUG
      void* old_val = p_ptr->exchange(b);
      BOOST_ASSERT (old_val == 0);
#else
      p_ptr->store(b);
#endif
    }
  };

  public:
  buffer_cache(transport trans, size_t buffer_size)
      : last_entry_written_plus_1(0), trans(trans), buffer_size(buffer_size), cache_deleted(new bool(false))
  {
    for (size_t i = 0; i < size; ++i) buffer_ptrs[i].store(0);
  }

  ~buffer_cache() {
    // The buffers themselves are owned by all_buffers so we don't need to free
    // them explicitly
    *cache_deleted = true;
  }

  boost::shared_ptr<void> allocate() {
    if (buffer_size == 0) return boost::shared_ptr<void>();
    void* p = 0;
    amplusplus::detail::atomic<void*>* p_ptr = 0;
    for (size_t i = 0; i < last_entry_written_plus_1.load(); ++i) {
      p = buffer_ptrs[i].exchange(0);
      if (p != 0) {
        p_ptr = &buffer_ptrs[i];
        goto have_buffer;
      }
    }
    {
      boost::shared_ptr<void> buf = trans.alloc_memory(buffer_size);
      BOOST_ASSERT (buf);
      {
        boost::lock_guard<amplusplus::detail::mutex> l(all_buffers_lock);
        all_buffers.push_back(buf);
        // fprintf(stderr, "%p increased buffer count to %zu\n", this, all_buffers.size());
      }
      p = buf.get();
      p_ptr = &buffer_ptrs[last_entry_written_plus_1.fetch_add(1)];
      if (p_ptr >= &buffer_ptrs[size]) {
        std::cerr << "Buffer list overflow" << std::endl;
        abort();
      }
    }
    have_buffer:
    BOOST_ASSERT (p);
    return boost::shared_ptr<void>(p, free_buf(cache_deleted, p_ptr));
  }
};
#endif

#if 0
class buffer_cache: boost::noncopyable {
  static const size_t size = (1 << 16);

  const size_t buffer_size;

  struct free_buf {
    void operator()(void* b) {
      AMPLUSPLUS_MPI_CALL_REGION MPI_Free_mem(b);
    }
  };

  public:
  buffer_cache(transport trans, size_t buffer_size)
      : buffer_size(buffer_size)
  {
  }

  ~buffer_cache() {}

  boost::shared_ptr<void> allocate() {
    if (buffer_size == 0) return boost::shared_ptr<void>();
    void* p = 0;
    AMPLUSPLUS_MPI_CALL_REGION MPI_Alloc_mem(buffer_size, MPI_INFO_NULL, &p);
    BOOST_ASSERT (p);
    return boost::shared_ptr<void>(p, free_buf());
  }
};
#endif

  }
}

#endif // AMPLUSPLUS_DETAIL_BUFFER_CACHE_HPP
