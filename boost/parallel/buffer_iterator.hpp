// Copyright (C) 2016-2017 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine 

#ifndef APPEND_BUFFER_ITERATOR_HPP
#define APPEND_BUFFER_ITERATOR_HPP

#include "append_buffer.hpp"

template <typename T, size_t LgFirstChunkSize = 21, size_t NumChunksToPrealloc = 1, size_t LgMaxSize = (sizeof(size_t) == 4 ? 31 : 40)>
class buffer_iterator: public std::iterator<std::random_access_iterator_tag, T> {

  typedef append_buffer<T, LgFirstChunkSize, NumChunksToPrealloc, LgMaxSize> Bucket;

public:
  buffer_iterator(): bucket(NULL), pos(0) {}
  buffer_iterator(Bucket* p, long po): bucket(p), pos(po) {}
  buffer_iterator(const buffer_iterator& other): bucket(other.bucket), pos(other.pos) {}

  buffer_iterator& operator++()    {pos++; return *this;} // prefix++
  buffer_iterator  operator++(int) {buffer_iterator tmp(*this); pos++; return tmp;} // postfix++
  buffer_iterator& operator--()    {
    //    assert(pos > 0);
    pos--; 
    return *this;
  } // prefix--

  buffer_iterator  operator--(int) {buffer_iterator tmp(*this); --(*this); return tmp;} // postfix--

  void     operator+=(const std::size_t& n)  {pos += n;}
  void     operator+=(const buffer_iterator& other) {pos += other.pos;}
  buffer_iterator operator+ (const std::size_t& n)  {buffer_iterator tmp(*this); tmp += n; return tmp;}
  buffer_iterator operator+ (const buffer_iterator& other) {buffer_iterator tmp(*this); tmp += other; return tmp;}

  void        operator-=(const std::size_t& n)  {
    //    assert(pos >= n);
    pos -= n;
  }
  void        operator-=(const buffer_iterator& other) {
    //    assert(pos >= other.pos);
    pos -= other.pos;
  }
  buffer_iterator    operator- (const std::size_t& n)  {buffer_iterator tmp(*this); tmp -= n; return tmp;}
  std::size_t operator- (const buffer_iterator& other) {
    //    assert(pos >= other.pos);
    return pos - other.pos;
  }

  bool operator< (const buffer_iterator& other) {return (pos-other.pos)< 0;}
  bool operator<=(const buffer_iterator& other) {return (pos-other.pos)<=0;}
  bool operator> (const buffer_iterator& other) {return (pos-other.pos)> 0;}
  bool operator>=(const buffer_iterator& other) {return (pos-other.pos)>=0;}
  bool operator==(const buffer_iterator& other) {return  pos == other.pos; }
  bool operator!=(const buffer_iterator& other) {return  pos != other.pos; }

  T& operator[](const int& n) {
    //fprintf(stderr, "============a\n");
    return (*bucket)[pos+n];
  }

  T& operator*() {
    //fprintf(stderr, "%zu============b\n", pos);
    return (*bucket)[pos];
  }

  T* operator->(){
    //fprintf(stderr, "============c\n");
    return &((*bucket)[pos]);
  }

  buffer_iterator operator=(const buffer_iterator& other) {
    pos = other.pos;	
    bucket = other.bucket; 
    return *this; 
  }

private:
  Bucket* bucket;
  long pos;

};

#endif
