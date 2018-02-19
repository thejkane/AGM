// Copyright (C) 2007-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/test/minimal.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/lexical_cast.hpp>
#include <ctime>

#define IS_MPI_TRANSPORT 1

using namespace boost;

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif


int test_main(int argc, char** argv)
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(1);

  int n = 10000;
  double p = 3e-3;
  int seed = std::time(0);

  if (argc > 1) n = lexical_cast<int>(argv[1]);
  if (argc > 2) p = lexical_cast<double>(argv[2]);
  if (argc > 3) seed = lexical_cast<int>(argv[3]);

  typedef adjacency_list<vecS, 
                         distributedS<vecS>,
                         directedS> Graph;

  int rank = trans.rank();
  int numprocs = trans.size();
  bool i_am_root = rank == 0;

  // broadcast the seed  
//   boost::mpi::broadcast(world, seed, 0);
  { amplusplus::scoped_epoch epoch(trans); }

  // Random number generator
  minstd_rand gen;

  if (i_am_root) {
    std::cout << "n = " << n << ", p = " << p << ", seed = " << seed
              << "\nBuilding graph with the iterator constructor... ";
    std::cout.flush();
  }

  parallel::variant_distribution<> distrib 
    = parallel::block<graph_traits<Graph>::vertices_size_type>(trans, n);

  // Build a graph using the iterator constructor, where each of the
  // processors has exactly the same information.
  gen.seed(seed);
  Graph g1(erdos_renyi_iterator<minstd_rand, Graph>(gen, n, p),
           erdos_renyi_iterator<minstd_rand, Graph>(),
           n, trans, distrib);

  { amplusplus::scoped_epoch epoch(trans); }

  sleep(5);

  // Build another, identical graph using add_edge
  if (i_am_root) {
    std::cout << "done.\nBuilding graph with add_edge from the root...";
    std::cout.flush();
  }

  gen.seed(seed);
  Graph g2(n, trans, distrib);

  { 
    amplusplus::scoped_epoch epoch(trans);

    if (i_am_root) {
      // The root will add all of the edges, some percentage of which
      // will require an immediate response from the owner of the edge.
      for (erdos_renyi_iterator<minstd_rand, Graph> first(gen, n, p), last;
	   first != last; ++first) {
	add_edge(vertex(first->first, g2), vertex(first->second, g2), g2);
	
	std::cout << rank << ": add_edge " << vertex(first->first, g2).local << "@" << vertex(first->first, g2).owner
		  << " -> " << vertex(first->second, g2).local << "@" << vertex(first->second, g2).owner << std::endl;
      }
    }
  }

  if (i_am_root) {
    std::cout << "synchronizing...";
    std::cout.flush();
  }

  // Verify that the two graphs are indeed identical.
  if (i_am_root) {
    std::cout << "done.\nVerifying graphs...";
    std::cout.flush();
  }

  // Check the number of vertices
  if (num_vertices(g1) != num_vertices(g2)) {
    std::cerr << g1.processor() << ": g1 has " << num_vertices(g1) 
              << " vertices, g2 has " << num_vertices(g2) << " vertices.\n";
    exit(-1);
  }

  // Check the number of edges
  if (num_edges(g1) != num_edges(g2)) {
    std::cerr << g1.processor() << ": g1 has " << num_edges(g1) 
              << " edges, g2 has " << num_edges(g2) << " edges.\n";
    exit(-1);
  }

  // Check the in-degree and out-degree of each vertex
  graph_traits<Graph>::vertex_iterator vfirst1, vlast1, vfirst2, vlast2;
  tie(vfirst1, vlast1) = vertices(g1);
  tie(vfirst2, vlast2) = vertices(g2);
  for(; vfirst1 != vlast1 && vfirst2 != vlast2; ++vfirst1, ++vfirst2) {
    if (out_degree(*vfirst1, g1) != out_degree(*vfirst2, g2)) {
      std::cerr << g1.processor() << ": out-degree mismatch ("
                << out_degree(*vfirst1, g1) << " vs. " 
                << out_degree(*vfirst2, g2) << ").\n";
      exit(-1);
    }
  }

  { amplusplus::scoped_epoch epoch(trans); }

  // TODO: Check the actual edge targets

  // Build another, identical graph using add_edge
  if (i_am_root) {
    std::cout << "done.\nBuilding graph with add_edge from everywhere...";
    std::cout.flush();
  }

  gen.seed(seed);
  Graph g3(n, trans, distrib);
  {
    // Each processor will take a chunk of incoming edges and add
    // them. Some percentage of the edge additions will require an
    // immediate response from the owner of the edge. This should
    // create a lot of traffic when building the graph, but should
    // produce a graph identical to the other graphs.
    typedef graph_traits<Graph>::edges_size_type edges_size_type;

    erdos_renyi_iterator<minstd_rand, Graph> first(gen, n, p);
    edges_size_type chunk_size = edges_size_type(p*n*n)/numprocs;
    edges_size_type start = chunk_size * rank;
    edges_size_type remaining_edges = 
      (rank < numprocs - 1? chunk_size 
       : edges_size_type(p*n*n) - start);

    // Spin the generator to the first edge we're responsible for
    for (; start; ++first, --start) 
      ;

    {
      amplusplus::scoped_epoch epoch(trans);

      for (; remaining_edges; --remaining_edges, ++first) {
	add_edge(vertex(first->first, g3), vertex(first->second, g3), g3);
      }
    }
  }

  if (i_am_root) {
    std::cout << "synchronizing...";
    std::cout.flush();
  }

  // Verify that the two graphs are indeed identical.
  if (i_am_root) {
    std::cout << "done.\nVerifying graphs...";
    std::cout.flush();
  }

  // Check the number of vertices
  if (num_vertices(g1) != num_vertices(g3)) {
    std::cerr << g1.processor() << ": g1 has " << num_vertices(g1) 
              << " vertices, g3 has " << num_vertices(g3) << " vertices.\n";
    exit(-1);
  }

  // Check the number of edges
  if (num_edges(g1) != num_edges(g3)) {
    std::cerr << g1.processor() << ": g1 has " << num_edges(g1) 
              << " edges, g3 has " << num_edges(g3) << " edges.\n";
    exit(-1);
  }

  // Check the in-degree and out-degree of each vertex
  tie(vfirst1, vlast1) = vertices(g1);
  tie(vfirst2, vlast2) = vertices(g3);
  for(; vfirst1 != vlast1 && vfirst2 != vlast2; ++vfirst1, ++vfirst2) {
    if (out_degree(*vfirst1, g1) != out_degree(*vfirst2, g3)) {
      std::cerr << g1.processor() << ": out-degree mismatch ("
                << out_degree(*vfirst1, g1) << " vs. " 
                << out_degree(*vfirst2, g3) << ").\n";
      exit(-1);
    }
  }

  // TODO: Check the actual edge targets

  if (i_am_root) {
    std::cout << "done.\n";
    std::cout.flush();
  }

  return 0;
}
