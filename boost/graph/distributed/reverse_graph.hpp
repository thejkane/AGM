// Copyright (C) 2005-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine
#ifndef BOOST_GRAPH_DISTRIBUTED_REVERSE_GRAPH_HPP
#define BOOST_GRAPH_DISTRIBUTED_REVERSE_GRAPH_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/parallel/container_traits.hpp>

namespace boost {
  namespace graph {

  /// Retrieve the process group from a reverse graph
  template<typename Graph, typename GraphRef>
  inline typename graph::parallel::process_group_type<Graph>::type
  transport(reverse_graph<Graph, GraphRef> const& g) {
    return g.m_g.transport();;

  }
} 

#endif
