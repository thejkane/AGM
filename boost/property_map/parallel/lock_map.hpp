// Copyright 2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine
#ifndef BOOST_PARALLEL_LOCK_MAP_HPP
#define BOOST_PARALLEL_LOCK_MAP_HPP

#include <boost/thread/mutex.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>

namespace boost { namespace parallel {

template <typename VertexIndexMap, typename Lock = boost::mutex>
class lock_map 
{
public:
  typedef typename property_traits<VertexIndexMap>::key_type     key_type;
  typedef typename property_traits<VertexIndexMap>::value_type   index_type;
  typedef boost::optional<shared_ptr<Lock> >                     value_type;
  typedef VertexIndexMap                                         vertex_index_map_type;
  typedef Lock                                                   lock_type;

  lock_map(VertexIndexMap vertex_index, index_type num_locks = 1)
      : vertex_index(vertex_index), locks(new Lock[num_locks]), num_locks(num_locks)
  { }

  // Condition below enables locking when false, intended for use with atomics_supported<key_type>
  template <typename Condition>
  typename boost::disable_if<Condition, value_type>::type
  maybe_lock(const key_type& key) const
  { return lock(key); }

  template <typename Condition>
  typename boost::enable_if<Condition, value_type>::type
  maybe_lock(const key_type&) const 
  { return value_type(); } // No lock returned

  // Lock unconditionally, return type always contains a lock
  value_type lock(const key_type& key) const
  { 
    locks[get(vertex_index, key) % num_locks].lock(); 
    return make_optional(true, boost::shared_ptr<Lock>(&locks[get(vertex_index, key) % num_locks], mem_fn(&Lock::unlock))); 
  }

private:
  VertexIndexMap            vertex_index;
  boost::shared_array<Lock> locks;
  index_type                num_locks;
};

template <typename Condition, typename VertexIndexMap, typename Lock>
typename lock_map<VertexIndexMap, Lock>::value_type
maybe_lock(lock_map<VertexIndexMap, Lock>& locks, typename lock_map<VertexIndexMap, Lock>::key_type& key)
{ return locks.template maybe_lock<Condition>(key); }

template <typename Condition, typename VertexIndexMap, typename Lock>
typename lock_map<VertexIndexMap, Lock>::value_type
maybe_lock(const lock_map<VertexIndexMap, Lock>& locks, typename lock_map<VertexIndexMap, Lock>::key_type& key)
{ return locks.template maybe_lock<Condition>(key); }

template <typename VertexIndexMap, typename Lock>
typename lock_map<VertexIndexMap, Lock>::value_type
lock(lock_map<VertexIndexMap, Lock>& locks, typename lock_map<VertexIndexMap, Lock>::key_type& key)
{ return locks.lock(key); }

template <typename VertexIndexMap, typename Lock = boost::mutex>
class dummy_lock_map 
{
public:
  typedef typename property_traits<VertexIndexMap>::key_type   key_type;
  typedef typename property_traits<VertexIndexMap>::value_type index_type;
  typedef boost::optional<boost::lock_guard<Lock> >            value_type;

  dummy_lock_map(VertexIndexMap) {}
  dummy_lock_map(VertexIndexMap, index_type) {}

  value_type maybe_lock(const key_type& key)
  { return value_type(); } // No lock returned

  value_type lock(const key_type& key)
  { return value_type(); } // No lock returned
};

template <typename VertexIndexMap, typename Lock>
typename dummy_lock_map<VertexIndexMap, Lock>::value_type
maybe_lock(dummy_lock_map<VertexIndexMap, Lock> locks, typename dummy_lock_map<VertexIndexMap, Lock>::key_type& key)
{ return locks.maybe_lock(key); }

template <typename VertexIndexMap, typename Lock>
typename dummy_lock_map<VertexIndexMap, Lock>::value_type
lock(dummy_lock_map<VertexIndexMap, Lock> locks, typename dummy_lock_map<VertexIndexMap, Lock>::key_type& key)
{ return locks.lock(key); }


} } // end namespace boost::parallel

#endif // BOOST_PARALLEL_LOCK_MAP_HPP
