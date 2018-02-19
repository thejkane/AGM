// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds

#include <am++/am++.hpp> 
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;

int main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  int num_threads = 2;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "--threads")
      num_threads = boost::lexical_cast<int>( argv[i+1] );
  }

  typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property, 
                                      distributedS<> >
    Digraph;

  // Build an Erdos-Renyi graph to test with
  typedef small_world_iterator<minstd_rand, Digraph> SWIter;

  int n = 50000, k = 4;
  double prob = 16. / n;

  minstd_rand gen;
  gen.seed(1);
  Digraph g(edges_are_unsorted, SWIter(gen, n, k, prob), SWIter(), n, trans);

  typedef two_bit_color_map<property_map<Digraph, vertex_index_t>::type> ColorMap;
  ColorMap color(num_vertices(g), get(vertex_index, g));

  boost::graph::distributed::breadth_first_search<Digraph> 
    bfs(g, make_bfs_visitor(null_visitor()), color);
  bfs.set_source(vertex(0, g));

  // Two threads now
  trans.set_nthreads(num_threads);

  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads]);
  for (int i = 0; i < num_threads; ++i) {
    boost::thread thr(boost::ref(bfs), i);
    threads[i].swap(thr);
  }

  for (int i = 0; i < num_threads; ++i)
    threads[i].join();

  // Back to one thread
  trans.set_nthreads(1);

  typedef property_traits<ColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;
  int visited = 0;
  BGL_FORALL_VERTICES(v, g, Digraph) { 
    if (get(color, v) == Color::black()) 
      ++visited; 
  }

  std::cerr << trans.rank() << ": visited " << visited << " vertices of " << num_vertices(g) << std::endl;

  return 0;
}
