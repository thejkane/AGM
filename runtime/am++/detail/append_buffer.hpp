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

#ifndef AMPLUSPLUS_DETAIL_APPEND_BUFFER_HPP
#define AMPLUSPLUS_DETAIL_APPEND_BUFFER_HPP

#include <boost/smart_ptr.hpp>
#include <boost/iterator.hpp>
#include <boost/assert.hpp>
#include <iterator>
#include <boost/intrusive/detail/utilities.hpp>
#include <cmath>
#include <am++/detail/thread_support.hpp>

// Container that can only grow in size, but has stable iterators and
// references.  It also has atomic push_back (i.e., multiple push_back
// operations do not conflict).  There is no guard that elements that have been
// pushed are actually stored before they are accessed (through iterators or
// operator[]()).  This container is a model of Random Access Container except
// for copy construction and assignment.  None of these operations are
// signal-safe (push_back(), operator[](), and iterator dereference can spin in
// some cases).

namespace amplusplus {
  namespace detail {

template <typename T>
class append_buffer;

template <typename T>
class append_buffer_iterator
    : public boost::iterator_facade<append_buffer_iterator<T>, T,
                                    std::random_access_iterator_tag>
{
  public:
  append_buffer_iterator(const append_buffer<T>& buf, size_t idx = 0)
    : ptr_to_current_chunk(&buf.chunk_ptrs[0]),
      current_chunk_num(0),
      offset_in_chunk(0),
      buf(buf)
  {this->advance((std::ptrdiff_t)idx);}

  private:
  T& dereference() const {
    T* ptr = 0;
    do {
      ptr = ptr_to_current_chunk->load();
    } while (!ptr); // Wait for chunk to be added if necessary
    return ptr[offset_in_chunk];
  }

  bool equal(const append_buffer_iterator<T>& o) const {
    return this->current_chunk_num == o.current_chunk_num &&
           this->offset_in_chunk == o.offset_in_chunk &&
           &this->buf == &o->buf;
  }

  void increment() {
    ++offset_in_chunk;
    if (offset_in_chunk == buf.chunk_size(current_chunk_num)) {
      ++current_chunk_num;
      ++ptr_to_current_chunk;
      offset_in_chunk = 0;
    }
  }

  void decrement() {
    if (offset_in_chunk == 0) {
      --current_chunk_num;
      --ptr_to_current_chunk;
      offset_in_chunk = buf.chunk_size(current_chunk_num);
    }
    --offset_in_chunk;
  }

  void advance(std::ptrdiff_t delta) {
    while (delta > 0) {
      size_t dist = (std::min)((size_t)delta, buf.chunk_size(current_chunk_num) - offset_in_chunk);
      offset_in_chunk += dist;
      delta -= dist;
      if (offset_in_chunk == buf.chunk_size(current_chunk_num)) {
        ++current_chunk_num;
        ++ptr_to_current_chunk;
        offset_in_chunk = 0;
      }
    }
    while (delta < 0) {
      if (offset_in_chunk == 0) {
        --current_chunk_num;
        --ptr_to_current_chunk;
        offset_in_chunk = buf.chunk_size(current_chunk_num);
      }
      size_t dist = (std::min)((size_t)(-delta), offset_in_chunk);
      offset_in_chunk -= dist;
      delta += dist;
    }
  }

  std::ptrdiff_t distance_to(const append_buffer_iterator<T>& o) const {
    return ((std::ptrdiff_t)(buf.chunk_offset(o.current_chunk_num)) -
              (std::ptrdiff_t)(buf.chunk_offset(this->current_chunk_num))) +
           ((std::ptrdiff_t)(o.offset_in_chunk) - (std::ptrdiff_t)(this->offset_in_chunk));
  }

  private:
  const amplusplus::detail::atomic<T*>* ptr_to_current_chunk;
  size_t current_chunk_num;
  size_t offset_in_chunk;
  const append_buffer<T>& buf;

  friend class iterator_core_access;
  friend class append_buffer<T>;
};

template <typename T>
class append_buffer {
  public:
  typedef T value_type;
  typedef append_buffer_iterator<T> iterator;
  typedef append_buffer_iterator<const T> const_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef std::ptrdiff_t difference_type;
  typedef size_t size_type;

  append_buffer(size_t initial_allocation = 16, size_t max_capacity = (size_t(1) << 31))
    : initial_allocation(initial_allocation),
      max_capacity(max_capacity),
      current_size(0),
      num_chunks(1 + ceil(log((double)max_capacity / initial_allocation) / log(2.))),
      chunk_ptrs(new amplusplus::detail::atomic<T*>[num_chunks])
  {
    // chunk_ptrs[0].store(new T[initial_allocation]);
    for(unsigned int i = 0; i < num_chunks; ++i)
      chunk_ptrs[i].store(0);
  }

  private: // Not implemented for now
  append_buffer(const append_buffer<T>&);
  append_buffer<T>& operator=(const append_buffer<T>&);

  public:
  ~append_buffer() {
    for (size_t i = 0; i < num_chunks; ++i) {
      delete[] chunk_ptrs[i].exchange(0);
    }
  }

  T& operator[](size_t idx) {
    size_t chunk_num = pos_to_chunk_num(idx);
    return chunk_ptrs[chunk_num][idx - chunk_offset(chunk_num)];
  }

  const T& operator[](size_t idx) const {
    size_t chunk_num = pos_to_chunk_num(idx);
    return chunk_ptrs[chunk_num][idx - chunk_offset(chunk_num)];
  }

  iterator begin() {
    return iterator(*this, 0);
  }

  const_iterator begin() const {
    return const_iterator(*this, 0);
  }

  iterator end() {
    return iterator(*this, current_size.load());
  }

  const_iterator end() const {
    return const_iterator(*this, current_size.load());
  }

  reverse_iterator rbegin() {
    return reverse_iterator(this->end());
  }

  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(this->end());
  }

  reverse_iterator rend() {
    return reverse_iterator(this->begin());
  }

  const_reverse_iterator rend() const {
    return reverse_iterator(this->begin());
  }

  size_t size() const {
    return current_size.load();
  }

  size_t max_size() const {
    return max_capacity;
  }

  bool empty() const {
    return size() == 0;
  }

  void swap(append_buffer<T>& o) {
    std::swap(initial_allocation, o.initial_allocation);
    std::swap(max_capacity, o.max_capacity);
    {
      const size_t temp = current_size.load();
      current_size.store(o.current_size.load());
      o.current_size.store(temp);
    }
    std::swap(num_chunks, o.num_chunks);
    chunk_ptrs.swap(o.chunk_ptrs);
  }

  bool operator==(const append_buffer<T>& o) const {
    return this->size() == o.size() &&
           std::equal(this->begin(), this->end(), o.begin());
  }

  bool operator!=(const append_buffer<T>& o) const {
    return !(*this == o);
  }

  bool operator<(const append_buffer<T>& o) const {
    return std::lexicographical_compare(this->begin(), this->end(), o.begin(), o.end());
  }

  bool operator>(const append_buffer<T>& o) const {
    return std::lexicographical_compare(o.begin(), o.end(), this->begin(), this->end());
  }

  bool operator<=(const append_buffer<T>& o) const {
    return !(*this > o);
  }

  bool operator>=(const append_buffer<T>& o) const {
    return !(*this < o);
  }

  size_t push_back_empty() { // Push back a default-constructed object
    size_t idx = current_size.fetch_add(1);
    BOOST_ASSERT (idx < max_capacity);
    size_t chunk_num = pos_to_chunk_num(idx);
    size_t offs = chunk_offset(chunk_num);
    size_t sz = chunk_size(chunk_num);
    amplusplus::detail::atomic<T*>& ptr = chunk_ptrs[chunk_num];
    if (idx == offs) { // I need to allocate this chunk
      T* this_chunk = new T[sz];
      T* old_chunk = ptr.exchange(this_chunk);
      BOOST_ASSERT (!old_chunk);
      return idx;
    } else { // Wait for the chunk to be allocated by someone else, then store into it
      T* chunk_start = 0;
      do {
        chunk_start = ptr.load();
      } while (!chunk_start);
      return idx;
    }
  }

  template <typename U>
  size_t push_back(const U& x) { // For anything that's assignable to T
    size_t loc = this->push_back_empty();
    (*this)[loc] = x;
    return loc;
  }

  private:
  size_t chunk_size(size_t idx) const {
    return initial_allocation << (idx == 0 ? 0 : idx - 1);
  }

  size_t chunk_offset(size_t idx) const {
    return (idx == 0) ? 0 : (initial_allocation << (idx - 1));
  }

  static size_t ffs_size_t(size_t x) {
    return (x == 0) ? 0 : 1 + boost::intrusive::detail::floor_log2(x);
  }

  size_t pos_to_chunk_num(size_t i) const {
    return ffs_size_t(i / initial_allocation);
  }

  private:
  size_t initial_allocation;
  size_t max_capacity;
  amplusplus::detail::atomic<size_t> current_size;
  size_t num_chunks;
  boost::shared_array<amplusplus::detail::atomic<T*> > chunk_ptrs;

  friend class append_buffer_iterator<T>;
};

  }
}

#endif // AMPLUSPLUS_DETAIL_APPEND_BUFFER_HPP
