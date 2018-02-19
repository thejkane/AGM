// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine
#ifndef BOOST_GRAPH_DISTRIBUTED_QUEUE_HPP
#define BOOST_GRAPH_DISTRIBUTED_QUEUE_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>

#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/parallel/graph_utility.hpp> // for identity
#include <algorithm> // for std::min, std::max

#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace boost { namespace graph { namespace distributed {

/// A unary predicate that always returns "true".                                                                                                                                                 
struct always_push
{
  template<typename T> bool operator()(const T&) const { return true; }
};

/** A distributed queue adaptor.
 *
 * Class template @c distributed_queue implements a distributed queue
 * across a process group. The distributed queue is an adaptor over an
 * existing (local) queue, which must model the @ref Buffer
 * concept. Each process stores a distinct copy of the local queue,
 * from which it draws or removes elements via the @ref pop and @ref
 * top members.
 *
 * The value type of the local queue must be a model of the @ref
 * GlobalDescriptor concept. The @ref push operation of the
 * distributed queue passes (via a message) the value to its owning
 * processor. Thus, the elements within a particular local queue are
 * guaranteed to have the process owning that local queue as an owner.
 *
 * Synchronization of distributed queues occurs in the @ref empty and
 * @ref size functions, which will only return "empty" values (true or
 * 0, respectively) when the entire distributed queue is empty. If the
 * local queue is empty but the distributed queue is not, the
 * operation will block until either condition changes. When the @ref
 * size function of a nonempty queue returns, it returns the size of
 * the local queue. These semantics were selected so that sequential
 * code that processes elements in the queue via the following idiom
 * can be parallelized via introduction of a distributed queue:
 *
 *   distributed_queue<...> Q;
 *   Q.push(x);
 *   while (!Q.empty()) {
 *     // do something, that may push a value onto Q
 *   }
 *
 * In the parallel version, the initial @ref push operation will place
 * the value @c x onto its owner's queue. All processes will
 * synchronize at the call to empty, and only the process owning @c x
 * will be allowed to execute the loop (@ref Q.empty() returns
 * false). This iteration may in turn push values onto other remote
 * queues, so when that process finishes execution of the loop body
 * and all processes synchronize again in @ref empty, more processes
 * may have nonempty local queues to execute. Once all local queues
 * are empty, @ref Q.empty() returns @c false for all processes.
 *
 * The distributed queue nearly models the @ref Buffer
 * concept. However, the @ref push routine does not necessarily
 * increase the result of @c size() by one (although the size of the
 * global queue does increase by one).
 */
size_t x;
template<typename OwnerMap, typename Buffer, 
         typename UnaryPredicate = always_push,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class distributed_queue
{
  typedef distributed_queue self_type;

 public:
  typedef Buffer                           buffer_type;
  typedef typename buffer_type::value_type value_type;
  typedef typename buffer_type::size_type  size_type;

  distributed_queue(amplusplus::transport& trans,
                    const OwnerMap& owner,
                    boost::shared_ptr<Buffer> buffer,
                    const UnaryPredicate& pred = UnaryPredicate(),
                    MessageGenerator message_gen = 
                      MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)));

  distributed_queue(amplusplus::transport& trans,
                    const OwnerMap& owner,
                    boost::shared_ptr<Buffer> buffer,
                    boost::shared_ptr<Buffer> incoming_buffer,
                    const UnaryPredicate& pred = UnaryPredicate(),
                    MessageGenerator message_gen = 
                      MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)));

  std::pair<unsigned long long, unsigned long long>
  get_cache_stats() {
    boost::parallel::all_reduce<unsigned long long, std::plus<unsigned long long> > 
      r(transport, std::plus<unsigned long long>());

    unsigned long long total_tests = r(push_msg.get()->counters.tests.load());
    unsigned long long total_hits = r(push_msg.get()->counters.hits.load());

    return std::make_pair(total_tests, total_hits);
  }

  /** Push an element onto the distributed queue.
   *
   * The element will be sent to its owner process to be added to that
   * process's local queue. 
   *
   * Complexity: O(1) messages of size O(sizeof(value_type)) will be
   * transmitted.
   */
  void push(const value_type& x) { push_msg.send(x); }

  /** Pop an element off the local queue. (NOT THREAD SAFE)
   *
   * @p @c !empty()
   */
  void pop() { buffer->pop(); }

  /**
   * Return the element at the top of the local queue. (NOT THREAD SAFE)
   *
   * @p @c !empty()
   */
  value_type top() { return buffer->top(); }

  /**
   * \overload (NOT THREAD SAFE)
   */
  const value_type& top() const { return buffer->top(); }

  /** Pop the top element off the queue and return it
   *
   */
  std::pair<value_type, bool> pop_top() { return buffer->pop_top(); }

  /** Pop the top 'N' elements off the queue and return them
   *
   *  std::vector is an arbitrary choice for the container returned
   */
  std::vector<value_type> pop_top_n(size_type n) { return buffer->pop_top_n(n); }

  /** Return a reference to the i'th element in the buffer
   *
   *  Combining this with a known size allows threaded iteration of
   *  all the elements in the queue in parallel with no atomics required
   *  for pop()
   */
  value_type& operator[](size_type i) { return (*buffer)[i]; }

  /** Clear queue (NOT THREAD SAFE)
   */
  void clear() { buffer->clear(); }

  /** Determine if the local queue is empty.
   *
   * TODO (NGE): Add all-reduce to determine when the whole
   *             distributed queue is empty.
   */
  bool empty() const { 
    using boost::parallel::all_reduce;

    all_reduce<size_type, std::plus<size_type> > r(transport, std::plus<size_type>());
    size_type total_size = r(buffer->size());
    return (total_size == 0); 
  }

  bool local_empty() const {return buffer->empty();}

  /** Determine the size of the local queue.
   */
  size_type size() const { return buffer->size(); }

  size_type incoming_size() const { return incoming_buffer->size(); } // DEBUG

  /** Swap current and incoming buffers if we're using a split-phase queue
   *
   *  (NOT THREAD SAFE)
   */
  void swap() { 
    assert(incoming_buffer); 
    boost::swap(buffer, incoming_buffer); 
  }

 private:

  // Message handlers
  struct push_handler;

  typedef typename MessageGenerator::template call_result<value_type, push_handler, OwnerMap,
				  amplusplus::duplicate_removal_t<boost::parallel::identity<value_type> > >::type
    push_message_type;

  const int dummy_first_member_for_init_order; // Unused
  amplusplus::transport& transport;
  OwnerMap owner;
  shared_ptr<Buffer> buffer;
  shared_ptr<Buffer> incoming_buffer;
  UnaryPredicate pred;

  push_message_type push_msg;
};

/// Helper macro containing the normal names for the template
/// parameters to distributed_queue.
#define BOOST_DISTRIBUTED_QUEUE_PARMS                           \
      typename OwnerMap, typename Buffer, typename UnaryPredicate, typename MessageGenerator

/// Helper macro containing the normal template-id for
/// distributed_queue.
#define BOOST_DISTRIBUTED_QUEUE_TYPE                                    \
      distributed_queue<OwnerMap, Buffer, UnaryPredicate, MessageGenerator>

/// Construct a new distributed queue.
template<typename OwnerMap, typename Buffer>
inline distributed_queue<OwnerMap, Buffer>
make_distributed_queue(amplusplus::transport& trans,
                       const OwnerMap& owner,
                       const Buffer& buffer)
{
  typedef distributed_queue<OwnerMap, Buffer> result_type;
  return result_type(trans, owner, buffer);
}

} } } // end namespace boost::graph::distributed

#include <boost/graph/distributed/detail/queue.cpp>

#undef BOOST_DISTRIBUTED_QUEUE_TYPE
#undef BOOST_DISTRIBUTED_QUEUE_PARMS

#endif // BOOST_GRAPH_DISTRIBUTED_QUEUE_HPP
