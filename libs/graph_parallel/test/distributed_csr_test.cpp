// Copyright (C) 2006-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

// A test of the distributed compressed sparse row graph type

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/test/minimal.hpp>

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

using namespace boost;

template <typename Transport>
void concept_checks(Transport& trans)
{
  typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property, 
                                      distributedS<> >
    Digraph;
  typedef graph_traits<Digraph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Digraph>::edge_descriptor edge_descriptor;

  function_requires< GraphConcept<Digraph> >();
  function_requires< IncidenceGraphConcept<Digraph> >();
  function_requires< AdjacencyGraphConcept<Digraph> >();

  function_requires< DistributedVertexListGraphConcept<Digraph> >();
  function_requires< DistributedEdgeListGraphConcept<Digraph> >();

  function_requires<
    ReadablePropertyGraphConcept<Digraph, vertex_descriptor, vertex_global_t> 
    >();
  function_requires<
    ReadablePropertyGraphConcept<Digraph, vertex_descriptor, vertex_owner_t>
    >();
  function_requires<
    ReadablePropertyGraphConcept<Digraph, vertex_descriptor, vertex_local_t>
    >();
  // NGE: There are some const issues in this concept check I believe
  function_requires<
    ReadablePropertyGraphConcept<Digraph, vertex_descriptor, vertex_index_t>
    >();

//   function_requires<
//     ReadablePropertyGraphConcept<Digraph, edge_descriptor, edge_global_t>
//     >();

  // DPG TBD: edge_owner, edge_local property maps

//   function_requires<
//     ReadablePropertyGraphConcept<Digraph, edge_descriptor, edge_index_t>
//     >();

  // Check default construction
  Digraph g(trans);
}

int test_main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(1);

  concept_checks(trans);

  typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property, 
                                      distributedS<> >
    Digraph;

  // Build an Erdos-Renyi graph to test with
  typedef sorted_erdos_renyi_iterator<minstd_rand, Digraph> ERIter;

  int n = 40;
  double prob = 0.1;

  minstd_rand gen;
  Digraph g(edges_are_unsorted, ERIter(gen, n, prob), ERIter(), n, trans);

  //  boost::graph::distributed::breadth_first_search<Digraph> bfs(g);
  typedef two_bit_color_map<property_map<Digraph, vertex_index_t>::type> ColorMap;
  ColorMap color(num_vertices(g), get(vertex_index, g));

  boost::graph::distributed::breadth_first_search<Digraph> 
    bfs(g, make_bfs_visitor(null_visitor()), color);

  bfs.run(vertex(0, g));

  typedef property_traits<ColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;
  BGL_FORALL_VERTICES(v, g, Digraph) { assert(get(color, v) == Color::black()); }

  std::ofstream out("dcsr.dot");
  write_graphviz(out, g);
  return 0;
}
