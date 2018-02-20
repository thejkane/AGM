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

#ifndef AMPLUSPLUS_DETAIL_VECTOR_OF_NONCOPYABLE_HPP
#define AMPLUSPLUS_DETAIL_VECTOR_OF_NONCOPYABLE_HPP

#include <boost/smart_ptr.hpp>
#include <boost/assert.hpp>
#include <utility>

namespace amplusplus {
  namespace detail {

template <typename T>
class vector_of_noncopyable {
  size_t size_, reserved_;
  boost::scoped_array<T> data_;

  void increase_size(size_t min_new_size) {
    size_t new_size = (std::max)(min_new_size, 2 * reserved_);
    boost::scoped_array<T> new_data(new T[new_size]);
    for (size_t i = 0; i < reserved_; ++i) {
      new_data[i].swap(data_[i]);
    }
    data_.swap(new_data);
    reserved_ = new_size;
  }

  public:
  vector_of_noncopyable(): data_(), size_(0), reserved_(0) {}
  vector_of_noncopyable(size_t sz): data_(new T[sz]), size_(sz), reserved_(sz) {}
  bool empty() const {return size_ == 0;}
  size_t size() const {return size_;}
  void push_back_empty() {
    if (size_ + 1 >= reserved_) increase_size(size_ + 1);
    BOOST_ASSERT (size_ < reserved_);
    ++size_;
  }
  void push_back_swap(T& t) {
    push_back_empty();
    data_[size_ - 1].swap(t);
  }
  T& back() {return data_.get()[size_ - 1];}
  const T& back() const {return data_.get()[size_ - 1];}
  T* begin() {return data_.get();}
  T* end() {return data_.get() + size_;}
  const T* begin() const {return data_.get();}
  const T* end() const {return data_.get() + size_;}
  T& operator[](size_t idx) {return data_[idx];}
  const T& operator[](size_t idx) const {return data_[idx];}
  void erase(T* iter) {
    --size_;
    T* end_ptr = end();
    for (; iter != end_ptr; ++iter) iter[0].swap(iter[1]);
    *end() = T();
    // FIXME: shrink
  }
};

  }
}

#endif // AMPLUSPLUS_DETAIL_VECTOR_OF_NONCOPYABLE_HPP
