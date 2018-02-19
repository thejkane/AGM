// Copyright 2011 The Trustees of Indiana University.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds

#ifndef BOOST_GRAPH_PARALLEL_ATOMIC_INSERT_QUEUE
#define BOOST_GRAPH_PARALLEL_ATOMIC_INSERT_QUEUE

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap

template <typename T>
class atomic_insert_queue
{
public:

  typedef T                  value_type;
  typedef unsigned long long size_type;

  atomic_insert_queue(size_type max_size)
    : buffer(new T[max_size]), max_size(max_size), head(0), tail(0) {}

  atomic_insert_queue(atomic_insert_queue& b)
    : buffer(new T[b.max_size]), max_size(b.max_size), head(0), tail(0) 
  { for (size_type i = 0; i < max_size; ++i) buffer[i] = b.buffer[i]; }

  atomic_insert_queue(const atomic_insert_queue& b)
    : buffer(new T[b.max_size]), max_size(b.max_size), head(0), tail(0)
  { for (size_type i = 0; i < max_size; ++i) buffer[i] = b.buffer[i]; }

  void push_back(const value_type& x) {
    using boost::parallel::fetch_and_add;
    buffer[fetch_and_add(&tail, (size_type)1) % max_size] = x; 
  }

  void pop() { 
    using boost::parallel::bool_compare_and_swap;

    if (head == tail) return;

    volatile size_type my_head = head;
    while (!bool_compare_and_swap(&head, my_head, my_head + 1)) {
      if (head >= tail) break;
      my_head = head;
    }
  }

  value_type& top() { return buffer[head % max_size]; }

  /** Pop and return the top element
   *
   *  The bool in the return value indicates whether the value_type is valid
   */
  std::pair<value_type, bool> pop_top() {
    using boost::parallel::bool_compare_and_swap;

    if (head < tail) {
      volatile size_type my_head = head;
      while (!bool_compare_and_swap(&head, my_head, my_head + 1)) {
	if (head >= tail) break;
	my_head = head;
      }
      return std::make_pair(buffer[my_head % max_size], true);
    }
    return std::make_pair(T(), false);
  }

  /** Return top 'N' elements
   *
   *  This may return 0-length vectors and is NOT deadlock-free
   */
  std::vector<value_type> pop_top_n(size_type n) { 

    // Getting the right number of elements here is the tricky part

    size_type s;
    while (head != tail) {
      using boost::parallel::bool_compare_and_swap;      

      volatile size_type my_head = head;
      s = std::min(n, (tail - head) % max_size);

      while (!bool_compare_and_swap(&head, my_head, my_head + s)) {
	if (head >= tail) break;
	my_head = head;
	s = std::min(n, (tail - head) % max_size);
      }
      return std::vector<value_type>(&buffer[my_head % max_size], 
				     &buffer[(my_head + s) % max_size]);
    }
    return std::vector<value_type>();
  }

  /** Return a reference to the i'th element in the buffer
   *
   *  Combining this with a known size allows threaded iteration of
   *  all the elements in the queue in parallel with no atomics required
   *  for pop()
   */
  value_type& operator[](size_type i) { return buffer[head + i % max_size]; }

  /** Determine if the queue is empty
   *
   *  This may return incorrect values, but nothing can actually be done 
   *  based on the return value anyway
   */
  bool empty() const { return (tail - head) == 0; }

  /** Determine the size of the local queue
   *
   *  This may return incorrect values, but nothing can actually be done 
   *  based on the return value anyway
   */
  size_type size() const { return tail - head; }

  size_type get_max_size() const { return max_size; }

  void clear() { head = tail = 0; }

protected:

  boost::scoped_array<T> buffer;
  size_type              max_size;
  size_type              head;
  size_type              tail;
};

#endif // BOOST_GRAPH_PARALLEL_ATOMIC_INSERT_QUEUE
