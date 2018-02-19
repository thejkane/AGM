// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>
#include <boost/graph/distributed/page_rank.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random.hpp>
#include <vector>
#include <iostream>
#include <stdlib.h>

#ifdef CSR
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#endif

using namespace boost;

bool close_to(double x, double y)
{
  double diff = x - y;
  if (diff < 0) diff = -diff;
  double base = (y == 0? x : y);
  if (base != 0) return diff / base < 0.01;
  else return true;
}

// Make convenient labels for the vertices
void test_distributed_page_rank(amplusplus::transport trans, int n, double p, int iterations, 
				int nthreads, bool verify, int seed)
{
  assert(nthreads == 1); // PageRank only works single threaded for now

#ifndef CSR
  typedef adjacency_list<vecS, 
                         distributedS<vecS>,
                         directedS> Graph;
#else
  typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property,
                                      distributedS<unsigned long> > Graph;
#endif

  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  minstd_rand gen;

  gen.seed(seed);

  Graph g(erdos_renyi_iterator<minstd_rand, Graph>(gen, n, p/2),
	  erdos_renyi_iterator<minstd_rand, Graph>(),
	  n, trans);

  // Setup RankMaps
  std::vector<double> rank(num_vertices(g)), other_rank(num_vertices(g));

  typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef property_map<Graph::base_type, vertex_index_t>::type LocalVertexIndexMap;
  typedef iterator_property_map<std::vector<double>::iterator, LocalVertexIndexMap> LocalRankMap;
  typedef property_map<Graph, vertex_global_t>::const_type VertexGlobalMap;
  typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap> RankMap;

  trans.set_nthreads(nthreads);
  boost::graph::distributed::page_rank<Graph, RankMap, boost::graph::n_iterations>
    PR(g, 
       make_iterator_property_map(rank.begin(), get(vertex_index, g)),
       make_iterator_property_map(other_rank.begin(), get(vertex_index, g)),
       boost::graph::n_iterations(iterations), 0.85, n);

  boost::scoped_array<boost::thread> threads(new boost::thread[nthreads]);
  for (int i = 0 ; i < nthreads ; ++i) {
    boost::thread thr(boost::ref(PR), i);
    threads[i].swap(thr);
  }
  
  for (int i = 0 ; i < nthreads ; ++i)
    threads[i].join();

  // Back to one thread
  trans.set_nthreads(1);

  if (trans.rank() == 0)
    std::cout << "Distributed PageRank complete (" << nthreads << " threads)\n";

  if ( verify ) {
    // Check against the sequential version

    typedef adjacency_list<listS, 
      vecS,
      directedS,
      // Vertex properties
      no_property,
      // Edge properties
      no_property > Graph2;
    
    gen.seed(seed);
    
    Graph2 g2(erdos_renyi_iterator<minstd_rand, Graph>(gen, n, p/2),
	      erdos_renyi_iterator<minstd_rand, Graph>(),
	      n);
    
    std::vector<double> rank2(num_vertices(g2));
    
    page_rank(g2, 
	      make_iterator_property_map(rank2.begin(), get(vertex_index, g2)),
	      boost::graph::n_iterations(iterations), 0.85, n);
    
    if (trans.rank() == 0)
      std::cout << "Sequential PageRank complete\n";

    VertexIndexMap v_index(get(vertex_index, g));
    
    boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
      global_index(trans, num_vertices(g), v_index, get(vertex_global, g));
    
    BGL_FORALL_VERTICES(v, g, Graph) {
      if (!close_to(rank[get(v_index, v)], rank2[get(global_index, v)]))
	std::cout << "distributed rank = " << rank[get(v_index, v)] << "  sequential rank = "
	          << rank2[get(global_index, v)] << std::endl;

      BOOST_CHECK(close_to(rank[get(v_index, v)], rank2[get(global_index, v)]));
    }
  }
}

int test_main(int argc, char* argv[])
{
  // AM++ initialization                                                                                                                                                                                                                       
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  int n = 1000;
  double p = 0.01;
  int iterations = 20;

  if (argc < 4) 
    test_distributed_page_rank(trans, n, p, iterations, 1, true, 123);
//     test_distributed_page_rank(trans, n, p, iterations, 2, true, 123);
  else 
    /* n, p, iterations threads seed verify */
    test_distributed_page_rank(trans, atoi(argv[1]), atof(argv[2]), atoi(argv[3]), 
			       argc == 4 ? 1 : atoi(argv[4]),
			       argc == 6 ? true : argv[6] == std::string("true"),
			       argc == 5 ? 123 : atoi(argv[5]));

  return 0;
}
