// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#ifndef BOOST_DISTRIBUTED_FILTERED_GRAPH_HPP
#define BOOST_DISTRIBUTED_FILTERED_GRAPH_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/graph/filtered_graph.hpp>

namespace boost {

  /// Retrieve the transport from a filtered graph
  template<typename Graph, typename EdgePredicate, typename VertexPredicate>
  inline amplusplus::transport&
  transport(filtered_graph<Graph, EdgePredicate, VertexPredicate> const& g) {
    return g.m_g.transport();
  }

  /// Forward vertex() to vertex() of the base graph 
  template <typename Graph, typename EdgePredicate, typename VertexPredicate>
  typename graph_traits<Graph>::vertex_descriptor
  vertex(typename graph_traits<Graph>::vertices_size_type i, 
         filtered_graph<Graph, EdgePredicate, VertexPredicate> const& g)
  { return vertex(i, g.m_g); }

  /// Forward make_vertex_descriptor() the base graph

}

#endif // BOOST_DISTRIBUTED_FILTERED_GRAPH_HPP
