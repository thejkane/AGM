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

#ifndef AMPLUSPLUS_LOCK_FREE_COALESCED_MESSAGE_TYPE_HPP
#define AMPLUSPLUS_LOCK_FREE_COALESCED_MESSAGE_TYPE_HPP

#include <am++/detail/atomics.hpp>
#include <am++/traits.hpp>
#include <utility>
#include <vector>
#include <boost/smart_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/assert.hpp>
#include <stdint.h>
#include <am++/scoped_raw_handler.hpp>

namespace amplusplus {

// Thread-safe functions:
//   Handler calls (from progress on underlying engine)
//   send
//   flush (both versions)
template <typename Arg, typename Handler, typename AMEngine, size_t BufferSize>
class lock_free_coalesced_message_type {
  struct message_buffer {
    volatile size_t elements_written;
    Arg data[BufferSize];
    message_buffer(): elements_written(0), data(data) {}
  };

  public:
  typedef Arg arg_type;
  typedef AMEngine engine_type;
  typedef Handler handler_type;
  typedef typename am_engine_traits<AMEngine>::rank_type rank_type;

  lock_free_coalesced_message_type(AMEngine& engine)
    : engine(engine),
      mpi_datatype(),
      handler_scope(engine, BufferSize, mpi_datatype.get(), BufferSize * sizeof(Arg), raw_message_handler(*this)),
      buffers_to_swap_in((void (*)(message_buffer*))0)
  {
    message_index = handler_scope.get_message_index();
    size_t nranks = num_processes(engine);
    for (size_t i = 0; i < 2 * nranks + 1; ++i) {
      buffer_pool.push_back(create_buffer());
    }
    for (size_t i = 0; i < nranks; ++i) {
      outgoing_buffers.push_back(std::make_pair(alloc_buffer(), 0));
    }
    flush_signal(engine, message_index).connect(typename lock_free_coalesced_message_type::flush_message_buffer(*this));
  }

  private:
  void send_buffer(message_buffer* buf,
                   uintptr_t count,
                   typename am_engine_traits<AMEngine>::rank_type dest) {
    while (buf->elements_written != count) {} // Wait for all elements to be written, not just allocated
    send(engine, message_index, buf->data, count, mpi_datatype.get(), dest, typename lock_free_coalesced_message_type::free_buffer(*this, buf));
  }

  struct flush_message_buffer {
    lock_free_coalesced_message_type& mt;
    flush_message_buffer(lock_free_coalesced_message_type& mt): mt(mt) {}
    void operator()(typename am_engine_traits<AMEngine>::rank_type dest) {
      BOOST_ASSERT (is_valid_rank(mt.engine, dest));
      if (mt.buffers_to_swap_in.get() == NULL) {
        mt.buffers_to_swap_in.reset(mt.alloc_buffer());
      }
      message_buffer* buffer_to_swap_in = mt.buffers_to_swap_in.get();
      volatile std::pair<message_buffer*, uintptr_t>* old_buf_ptr = &mt.outgoing_buffers[dest];
      message_buffer* old_buf_first = 0;
      uintptr_t old_buf_second = 0;
      detail::atomic_dwcas(old_buf_ptr, (message_buffer*)0, (uintptr_t)0, (message_buffer*)0, (uintptr_t)0, old_buf_first, old_buf_second); // Read atomically
      while (true) { // Use CAS for atomic swap
        if (old_buf_second == 0) return;
        message_buffer* new_buf_first = buffer_to_swap_in;
        uintptr_t new_buf_second = 0;
        bool worked = detail::atomic_dwcas(old_buf_ptr, old_buf_first, old_buf_second, new_buf_first, new_buf_second, old_buf_first, old_buf_second);
        if (worked) break;
      }
      uintptr_t count = old_buf_second;
      mt.send_buffer(old_buf_first, count, dest);
      mt.buffers_to_swap_in.reset(mt.alloc_buffer());
    }
  };

  struct raw_message_handler {
    lock_free_coalesced_message_type& mt;
    raw_message_handler(lock_free_coalesced_message_type& mt): mt(mt) {}
    void operator()(typename am_engine_traits<AMEngine>::rank_type src, const void* buf, size_t count) {
      const Arg* actual_buf = (const Arg*)buf;
      for (size_t i = 0; i < count; ++i) {
        mt.handler(src, actual_buf[i]);
      }
    }
  };
  
  public:
  friend void set_handler(lock_free_coalesced_message_type& mt, const handler_type& handler) {
    mt.handler = handler;
  }

  friend void send(lock_free_coalesced_message_type& mt, const arg_type& arg, rank_type dest) {
    BOOST_ASSERT (is_valid_rank(mt.engine, dest));
    if (mt.buffers_to_swap_in.get() == NULL) {
      mt.buffers_to_swap_in.reset(mt.alloc_buffer());
    }
    message_buffer* buffer_to_swap_in = mt.buffers_to_swap_in.get();
    volatile std::pair<message_buffer*, uintptr_t>* old_buf_ptr = &mt.outgoing_buffers[dest];
    message_buffer* old_buf_first = 0;
    uintptr_t old_buf_second = 0;
    detail::atomic_dwcas(old_buf_ptr, (message_buffer*)0, (uintptr_t)0, (message_buffer*)0, (uintptr_t)0, old_buf_first, old_buf_second); // Read atomically
    while (true) {
      message_buffer* new_buf_first = old_buf_second + 1 == BufferSize ? buffer_to_swap_in : old_buf_first;
      uintptr_t new_buf_second = old_buf_second + 1 == BufferSize ? 0 : old_buf_second + 1;
      bool worked = detail::atomic_dwcas(old_buf_ptr, old_buf_first, old_buf_second, new_buf_first, new_buf_second, old_buf_first, old_buf_second);
      if (worked) break;
    }
    uintptr_t count = old_buf_second;
    old_buf_first->data[count] = arg;
    detail::atomic_inc(old_buf_first->elements_written);
    if (count + 1 == BufferSize) {
      mt.send_buffer(old_buf_first, count + 1, dest);
      mt.buffers_to_swap_in.reset(mt.alloc_buffer());
    }
  }

  friend void flush(lock_free_coalesced_message_type& mt) {
    flush(mt.engine, mt.message_index); // Calls my flush through hook signal
  }

  friend void flush(lock_free_coalesced_message_type& mt, typename am_engine_traits<AMEngine>::rank_type dest) {
    flush(mt.engine, dest, mt.message_index); // Calls my flush through hook signal
  }

  private:
  AMEngine& engine;
  make_mpi_datatype<Arg> mpi_datatype;
  typename am_engine_traits<AMEngine>::msg_index_type message_index;
  handler_type handler;
  scoped_raw_handler<AMEngine> handler_scope;
  std::vector<boost::shared_ptr<void> > all_buffers; // Free and used, keeps ownership of them
  std::vector<message_buffer*> buffer_pool;
  std::vector<std::pair<message_buffer*, uintptr_t> > outgoing_buffers; // Pair is <buffer, count>
  boost::thread_specific_ptr<message_buffer> buffers_to_swap_in; // One per thread

  message_buffer* create_buffer() {
    boost::shared_ptr<void> buf = alloc_memory(engine, sizeof(message_buffer));
    all_buffers.push_back(buf);
    message_buffer* result = (message_buffer*)buf.get();
    result->elements_written = 0;
    return result;
  }

  message_buffer* alloc_buffer() {
    if (buffer_pool.empty()) {
      return create_buffer();
    } else {
      message_buffer* buf = buffer_pool.back();
      buffer_pool.pop_back();
      buf->elements_written = 0;
      return buf;
    }
  }

  struct free_buffer {
    lock_free_coalesced_message_type& mt;
    message_buffer* buf;
    free_buffer(lock_free_coalesced_message_type& mt, message_buffer* buf): mt(mt), buf(buf) {}
    void operator()(const MPI_Status&) const {
      mt.buffer_pool.push_back(buf);
    }
  };
};

}

#endif // AMPLUSPLUS_LOCK_FREE_COALESCED_MESSAGE_TYPE_HPP
