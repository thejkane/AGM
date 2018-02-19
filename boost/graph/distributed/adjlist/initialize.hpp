// Copyright (C) 2007-2012 The Trustees of Indiana University.
// Copyright (C) 2007 Douglas Gregor 

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Nicholas Edmonds

// This file contains code for the distributed adjacency list's
// initializations. It should not be included directly by users.

#ifndef BOOST_GRAPH_DISTRIBUTED_ADJLIST_INITIALIZE_HPP
#define BOOST_GRAPH_DISTRIBUTED_ADJLIST_INITIALIZE_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

namespace boost { 

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
template<typename EdgeIterator>
void
PBGL_DISTRIB_ADJLIST_TYPE::
initialize(EdgeIterator first, EdgeIterator last,
           vertices_size_type, const base_distribution_type& distribution, 
           vecS)
{
  rank_type id = transport_.rank();
  while (first != last) {
    if ((rank_type)distribution(first->first) == id) {
      vertex_descriptor source(id, distribution.local(first->first));
      vertex_descriptor target(distribution(first->second),
                               distribution.local(first->second));
      add_edge(source, target, *this);
    }
    ++first;
  }
}

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
template<typename EdgeIterator, typename EdgePropertyIterator>
void
PBGL_DISTRIB_ADJLIST_TYPE::
initialize(EdgeIterator first, EdgeIterator last,
           EdgePropertyIterator ep_iter,
           vertices_size_type, const base_distribution_type& distribution, 
           vecS)
{
  rank_type id = transport_.rank();
  while (first != last) {
    if (static_cast<rank_type>(distribution(first->first)) == id) {
      vertex_descriptor source(id, distribution.local(first->first));
      vertex_descriptor target(distribution(first->second),
                               distribution.local(first->second));
      add_edge(source, target, *ep_iter, *this);
    }
    ++first;
    ++ep_iter;
  }
}

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
template<typename EdgeIterator, typename EdgePropertyIterator,
         typename VertexListS>
void
PBGL_DISTRIB_ADJLIST_TYPE::
initialize(EdgeIterator first, EdgeIterator last,
           EdgePropertyIterator ep_iter,
           vertices_size_type n, const base_distribution_type& distribution,
           VertexListS)
{
  assert(false); // NGE: ListS for vertices doesn't work at the moment
}

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
template<typename EdgeIterator, typename VertexListS>
void
PBGL_DISTRIB_ADJLIST_TYPE::
initialize(EdgeIterator first, EdgeIterator last,
           vertices_size_type n, const base_distribution_type& distribution,
           VertexListS)
{
  assert(false); // NGE: ListS for vertices doesn't work at the moment
}

}   // end namespace boost

#endif // BOOST_GRAPH_DISTRIBUTED_ADJLIST_INITIALIZE_HPP
