// Copyright (C) 2007-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GRAPH_OWNER_DEFS
#define BOOST_GRAPH_OWNER_DEFS

#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/format.hpp>
#include "boost/tuple/tuple.hpp"
//#include <boost/tuple/tuple_comparison.hpp>

namespace {
  template <typename OwnerMap, typename Data>
  struct owner_from_tuple {
    explicit owner_from_tuple(const OwnerMap& owner) : owner(owner) {}
    
    const OwnerMap& owner;
  };
  
  template <typename OwnerMap, typename Data>
  typename boost::property_traits<OwnerMap>::value_type
  get(const owner_from_tuple<OwnerMap, Data>& o, const Data& data)
  { return get(o.owner, data.template get<0>()); }

  template <typename OwnerMap, typename Data>
  struct owner_from_pair {
    explicit owner_from_pair(const OwnerMap& owner) : owner(owner) {}
    
    const OwnerMap& owner;
  };
  
  template <typename OwnerMap, typename Data>
  typename boost::property_traits<OwnerMap>::value_type
  get(const owner_from_pair<OwnerMap, Data>& o, const Data& data)
  { return get(o.owner, data.first); }

  template <typename OwnerMap, typename PredecessorMap, typename Data>
  struct owner_from_predecessor {
    explicit owner_from_predecessor(const OwnerMap& owner, const PredecessorMap& predecessor) : owner(owner), predecessor(predecessor) {}
    
    const OwnerMap& owner;
    const PredecessorMap& predecessor;
  };

  template <typename OwnerMap, typename PredecessorMap, typename Data>
  typename boost::property_traits<OwnerMap>::value_type
  get(const owner_from_predecessor<OwnerMap, PredecessorMap, Data>& o, const Data& data) 
  { return get(o.owner, get(o.predecessor, data.template get<0>())); }

  struct minimum_tuple_first
  {
    template<typename T>
    const T& operator()(const T& x, const T& y) const { return x.template get<0>() < y.template get<0>() ? x : y; }

    template<typename F>
    struct result {
      typedef typename boost::function_traits<F>::arg1_type type;
    };
  };
}

#endif
