// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine 

#ifndef APPEND_BUFFER_HPP
#define APPEND_BUFFER_HPP

#ifndef __GNU_SOURCE
#define __GNU_SOURCE
#endif

#include <vector>
#include <utility>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <am++/detail/thread_support.hpp>
#include <string.h>
#include <limits.h>

// Buffer with append and block access capabilities.  The buffer is meant to be
// used in the following pattern:
//
// 1. Create buffer.
// 2. Add elements (using append()).
// 3. Access elements (using get_data() and size()).
// 4. [Optionally] Clear buffer using clear() and go to step 2.
//
// T is the element type; it should be default and copy constructible.
// LgFirstChunkSize is the logarithm base 2 of the size of the first (smallest)
// chunk of data; LgMaxSize is the log base 2 of the maximum size of the
// buffer.  The buffer was not tested with LgFirstChunkSize == 0 or LgMaxSize
// near the number of bits in a size_t; it is unlikely to work in those cases.
// In particular, some intrinsics assume signed numbers while size_t is
// unsigned, so doing anything that would use the sign bit of a size may not
// work.
//
// To avoid data copying, the buffer is arranged as an array of chunks; the
// chunks get exponentially larger as more elements are inserted.  The first
// and second chunks have size 2**LgFirstChunkSize, and the sizes double from
// there.  The get_data() method returns a vector of thunk base pointers and
// offsets; the user is not responsible to manage them (clear() and the
// destructor do that).  Even though two loops are required to access the
// elements, the smallest inner loop has size 2**LgFirstChunkSize and they
// become much larger if there are many elements in the queue.  Note that
// append() is marked with parallel loops but with a count of 1; they might
// need to be marked sequential if the compiler tries to parallelize them
// despite the count.  The common case is for a write to be within a particular
// chunk that does not need to be allocated; little of the code in append() is
// used in that case.
//
template <typename T, size_t LgFirstChunkSize = 21, size_t NumChunksToPrealloc = 1, size_t LgMaxSize = (sizeof(size_t) == 4 ? 31 : 40)>
class append_buffer {
  public:

  typedef T value_type;
  typedef size_t size_type;

  append_buffer() {
    current_size.store(0);
    for (size_t i = 0; i < NumChunksToPrealloc; ++i) {
      chunks[i].store(new T[chunk_size(i)]);
    }
    for (size_t i = NumChunksToPrealloc; i < NumChunks; ++i) {
      chunks[i].store(0);
    }
  }

  size_t size() const {return current_size.load();}
  bool empty() const {return current_size.load() == 0;}

  template <typename Iter>
  void append(Iter data, size_t count) {
    if (count == 0) return;
    size_t first_index = current_size.fetch_add(count);
    size_t last_index = first_index + count;
    // fprintf(stderr, "Copying into logical %zu..%zu\n", first_index, last_index);
    // Do any allocations that need to be done in this append (the first write
    // to a chunk must allocate it).  Do this first because other writes may be
    // blocked waiting for this allocation.
    size_t first_chunk_written = chunk_for_index(first_index);
    size_t last_chunk_written = chunk_for_index(last_index - 1);
    // fprintf(stderr, "Writing from chunks %zu to %zu\n", first_chunk_written, last_chunk_written);

    for (size_t i = (std::max)(NumChunksToPrealloc, first_chunk_written); i <= last_chunk_written; ++i) {
      size_t this_chunk_start = chunk_first_index(i);
      if (this_chunk_start >= first_index && this_chunk_start < last_index) {
        // fprintf(stderr, "Allocating chunk %zu of size %zu\n", i, chunk_size(i));
        T* this_chunk = new T[chunk_size(i)];
#ifdef NDEBUG
        chunks[i].store(this_chunk);
#else
        T* old_val = chunks[i].exchange(this_chunk);
        assert (old_val == 0);
#endif
      }
    }
    // Do the writes

    for (size_t i = first_chunk_written; i <= last_chunk_written; ++i) {
      size_t this_chunk_start = chunk_first_index(i);
      size_t this_chunk_size = chunk_size(i);
      size_t this_offset_in_data = (std::max)(first_index, this_chunk_start) - first_index;
      size_t this_offset_in_chunk = (std::max)(first_index, this_chunk_start) - this_chunk_start;
      size_t this_size_to_write = (std::min)(this_chunk_size - this_offset_in_chunk, count - this_offset_in_data);
      // fprintf(stderr, "Copying from offset %zu in data buffer to offset %zu in chunk %zu; size = %zu\n", this_offset_in_data, this_offset_in_chunk, i, this_size_to_write);
      const T* __restrict__ this_source_ptr = data + this_offset_in_data;
      // This statement may block if the chunk has not been allocated yet.
      T* start_ptr = chunks[i].load();
      while (start_ptr == 0) {start_ptr = chunks[i].load();} // Wait for chunk to be allocated if necessary
      T* __restrict__ this_dest_ptr = start_ptr + this_offset_in_chunk;
      std::copy(this_source_ptr, this_source_ptr + this_size_to_write, this_dest_ptr);
    }
  }

  T& operator[](size_t index) {
    size_t chunk_num = chunk_for_index(index);
    T* chunk_ptr = chunks[chunk_num].load();
    size_t chunk_start = chunk_first_index(chunk_num);
    T* __restrict__ value_ptr = chunk_ptr + index - chunk_start;
    return *value_ptr;
  }

  void push_back(T data) {
    size_t first_index = current_size.fetch_add(1);
    // fprintf(stderr, "Copying into logical %zu..%zu\n", first_index, last_index);
    // Do any allocations that need to be done in this append (the first write
    // to a chunk must allocate it).  Do this first because other writes may be
    // blocked waiting for this allocation.
    size_t chunk_written = chunk_for_index(first_index);
    size_t this_chunk_start = chunk_first_index(chunk_written);
    if (this_chunk_start == first_index && chunk_written >= NumChunksToPrealloc) {
      // fprintf(stderr, "Allocating chunk %zu of size %zu\n", chunk_written, chunk_size(chunk_written));
      T* this_chunk = new T[chunk_size(chunk_written)];
#ifdef NDEBUG
      chunks[chunk_written].store(this_chunk);
#else
      T* old_val = chunks[chunk_written].exchange(this_chunk);
      assert (old_val == 0);
#endif
    }
    // Do the write
    // This statement may block if the chunk has not been allocated yet.
    T* start_ptr = chunks[chunk_written].load();
    while (start_ptr == 0) {start_ptr = chunks[chunk_written].load();} // Wait for chunk to be allocated if necessary
    T* __restrict__ this_dest_ptr = start_ptr + first_index - this_chunk_start;
    *this_dest_ptr = data;
  }

  template <size_t LgFirstChunkSize2, size_t LgMaxSize2>
  void append(const append_buffer<T, LgFirstChunkSize2, LgMaxSize2>& other_buf) {
    std::vector<std::pair<const T* __restrict__, size_t> > ranges;
    other_buf.get_data(ranges);
    size_t rsize = ranges.size();
    for (size_t i = 0; i < rsize; ++i) {
      append(ranges[i].first, ranges[i].second);
    }
  }

  template <typename Iter>
  void push_back(Iter a, Iter b) {
    append(a, std::distance(a, b));
  }

  void get_data(std::vector<std::pair<const T* __restrict__, size_t> >& result) const {
    size_t old_result_size = result.size();
    size_t num_chunks = num_chunks_in_use();
    size_t cur_size = current_size;
    result.resize(num_chunks + old_result_size);
    for (size_t i = 0; i < num_chunks; ++i) {
      result[old_result_size + i].first = chunks[i].load();
      result[old_result_size + i].second =
        (std::min)(chunk_size(i), cur_size - chunk_first_index(i));
    }
  }

  void clear() {
    size_t num_chunks = num_chunks_in_use();

    for (size_t i = NumChunksToPrealloc; i < num_chunks; ++i) {
      delete[] chunks[i].exchange(0);
    }
    current_size.store(0);
  }

  void swap(append_buffer& b) {
    std::swap(current_size, b.current_size);
    for (size_t i = 0; i < NumChunks; ++i) {T* p = chunks[i].load(); chunks[i].store(b.chunks[i].exchange(p));}
  }

  ~append_buffer() {
    for (size_t i = 0; i < NumChunks; ++i) {
      delete[] chunks[i].exchange(0);
    }
    current_size.store(0);
  }

  private: // Helper definitions and functions
  static const size_t FirstChunkSize = ((size_t)1) << LgFirstChunkSize;
  static const size_t MaxSize = ((size_t)1) << LgMaxSize;
  static const size_t NumChunks = LgMaxSize - LgFirstChunkSize + 1;

  static inline size_t chunk_size(size_t i) {
    assert (i < NumChunks);
    return ((size_t)1) << ((std::max)(i, (size_t)1) + LgFirstChunkSize - 1);
  }

  static inline size_t chunk_first_index(size_t i) {
    assert (i < NumChunks);
    return i == 0 ? 0 : (((size_t)1) << (i + LgFirstChunkSize - 1));
  }

  static inline size_t chunk_for_index(size_t i) {
    assert (i < MaxSize);
    size_t bits_used = (i == 0) ? 0 : (sizeof(size_t) * CHAR_BIT - __builtin_clzl(i));
    size_t result = (std::max)(bits_used, LgFirstChunkSize) - LgFirstChunkSize;
    assert (result < NumChunks);
    assert (chunk_first_index(result + 1) == chunk_first_index(result) + chunk_size(result));
    assert (i >= chunk_first_index(result) && i < chunk_first_index(result + 1));
    // fprintf(stderr, "chunk_for_index(%zu) = %zu\n", i, result);
    return result;
  }

  inline size_t num_chunks_in_use() const {
    size_t cur_size = current_size.load();
    return cur_size == 0 ? 0 : chunk_for_index(cur_size - 1) + 1;
  }

  private: // Data members
  // Current size of buffer; always kept full for use with int_fetch_add.  Must
  // be updated atomically.
  amplusplus::detail::atomic<size_t> current_size;
  // Chunks of data.  Unallocated chunk pointers are NULL and empty so trying
  // to access unallocated chunks will block.  The append() method determines
  // itself which thread needs to allocate a particular chunk.
  amplusplus::detail::atomic<T*> chunks[NumChunks];
};

#endif // APPEND_BUFFER_HPP
