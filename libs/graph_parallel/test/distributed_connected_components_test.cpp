// Copyright (C) 2004-2006 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#define PRINT_DEBUG

#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/connected_components.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/distributed/graphviz.hpp>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <boost/random.hpp>
#include <boost/test/minimal.hpp>

#ifdef CSR
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#endif

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

using namespace boost;

typedef double time_type;

inline time_type get_time()
{
        return MPI_Wtime();
}

std::string print_time(time_type t)
{
  std::ostringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
  return out.str();
}

template<typename T>
class map_lt
{
public:
  bool operator()() const { return false; }
  bool operator()(T x, T y) const { return (owner(x) < owner(y) || (owner(x) == owner(y) && local(x) < local(y))); }
};

void 
test_distributed_connected_components(amplusplus::transport trans, int n, double _p, bool verify, 
				      bool emit_dot_file, int seed, int nthreads)
{
#ifndef CSR
  typedef adjacency_list<vecS, 
                         distributedS<vecS>,
                         directedS> Graph;
#else
  typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property,
                                      distributedS<unsigned long> > Graph;
#endif

  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits<Graph>::vertices_size_type vertices_size_type;

  typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef vector_property_map<Vertex, VertexIndexMap> ParentMap;
  typedef property_map<Graph, vertex_global_t>::const_type VertexGlobalMap;

  minstd_rand gen;

  gen.seed(seed);

  parallel::variant_distribution<vertices_size_type> distrib 
    = parallel::oned_block_cyclic<vertices_size_type>(trans, 2);

  typedef graph_traits<Graph>::edges_size_type edges_size_type;
  std::vector<std::pair<edges_size_type, edges_size_type> > edges;

  typedef erdos_renyi_iterator<minstd_rand, Graph> ERIter;

  for (ERIter edges_begin(gen, n, _p) ; edges_begin != ERIter() ; ++edges_begin) {
    if (distrib(edges_begin->first) == trans.rank())
      edges.push_back(std::make_pair(edges_begin->first, edges_begin->second));
    if (distrib(edges_begin->second) == trans.rank())
      edges.push_back(std::make_pair(edges_begin->second, edges_begin->first));
  }

#ifdef CSR
  Graph g(edges_are_unsorted, edges.begin(), edges.end(), n, trans, distrib);
#else
  // Use directedS and flip edges explicitly for now
  Graph g(edges.begin(), edges.end(), n, trans, distrib);
#endif

  edges.clear();

  // ParentMap
  ParentMap parent(num_vertices(g), get(vertex_index, g));

  BGL_FORALL_VERTICES(v, g, Graph) {
    put(parent, v, v); // Initialize to self because we're going to use null_vertex() as a sentinel value in parallel search
  }

  // DEBUG
  typedef boost::parallel::lock_map<VertexIndexMap> LockMap;
  LockMap locks(get(vertex_index, g), num_vertices(g) / 64);

  typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
  MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(4096)));
  // DEBUG

  // nthreads threads in transport now
  trans.set_nthreads(nthreads);
  boost::graph::distributed::connected_components<Graph, ParentMap, MessageGenerator, LockMap> CC(g, parent, locks);
	
  boost::scoped_array<boost::thread> threads(new boost::thread[nthreads]);
  for (int i = 0; i < nthreads; ++i) {
    boost::thread thr(boost::ref(CC), i);
    threads[i].swap(thr);
  }

  for (int i = 0; i < nthreads; ++i)
    threads[i].join();

  // Back to one thread
  trans.set_nthreads(1);

  BGL_FORALL_VERTICES(v, g, Graph) {
    BGL_FORALL_ADJ(v, u, g, Graph) {
      if (get(get(vertex_owner, g), u) == trans.rank())
	assert(get(parent, v) == get(parent, u));
    }
  }

  // Number components
  typedef property_map<Graph::base_type, vertex_index_t>::type LocalVertexIndexMap;
  typedef vector_property_map<int, LocalVertexIndexMap> LocalComponentMap;
  typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalComponentMap> ComponentMap;
  ComponentMap component(trans, get(vertex_global, g), LocalComponentMap(num_vertices(g), get(vertex_index, g.base())));

  int num_components = CC.number_components<ComponentMap>(component);

  if (trans.rank() == 0)
    std::cout << "Distributed CC found " << num_components << " components\n";

  if ( verify )
    {
      if ( trans.rank() == 0 )
        {
          component.set_max_ghost_cells(0);

	  {
	    amplusplus::scoped_epoch epoch(trans);

	    BGL_FORALL_VERTICES(v, g, Graph) {
	      get(component, v);
	    }
	  }

          // Check against the sequential version
          typedef adjacency_list<listS, 
            vecS,
            undirectedS,
            // Vertex properties
            no_property,
            // Edge properties
            no_property > Graph2;
          
          gen.seed(seed);
          
	  Graph2 g2(erdos_renyi_iterator<minstd_rand, Graph>(gen, n, _p/2),
		    erdos_renyi_iterator<minstd_rand, Graph>(),
		    n);

          std::vector<int> component2 (n);
          int tmp;
          tmp = connected_components(g2, make_iterator_property_map(component2.begin(), get(vertex_index, g2)));
          std::cout << "Verifier found " << tmp << " components" << std::endl;
          
          // Make sure components and component2 match
          std::map<int, int> c2c;
          int i;
          // This fails if there are more components in 'component' than
          // 'component2' because multiple components in 'component' may 
          // legitimately map to the same component number in 'component2'.
          // We can either test the mapping in the opposite direction or
          // just assert that the numbers of components found by both 
          // algorithms is the same
	  {
	    amplusplus::scoped_epoch epoch(trans); // distributed_property_map::get() sends messages

	    for ( i = 0; i < n; i++ )
	      if ( c2c.find( get(component, vertex(i, g)) ) == c2c.end() )
		c2c[get(component, vertex(i, g))] = component2[i];
	      else
		if ( c2c[get(component, vertex(i, g))] != component2[i] )
		  break;
	  }
	  
	  if ( i < n || num_components != tmp) {
	    printf("Unable to verify CC result...\n");
          } else
            printf("Passed verification... %i connected components\n", 
                   (int)c2c.size());
        }
      else 
        {
	  { amplusplus::scoped_epoch epoch(trans); }
	  { amplusplus::scoped_epoch epoch(trans); }
        }
      if ( emit_dot_file )
        write_graphviz("cc.dot", g, paint_by_number(component));
    }
}

int test_main(int argc, char* argv[])
{
  // AM++ initialization                                                                                                                                                                                                                       
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  if ( argc < 6 ) {
    test_distributed_connected_components(trans, 10000, 0.001, true, false, 1, 1);
    test_distributed_connected_components(trans, 10000, 0.001, true, false, 1, 2);
  }
  else
    test_distributed_connected_components
      (trans, atoi(argv[1]), atof(argv[2]), 
       argv[3]==std::string("true"), argv[4]==std::string("true"),
       argc == 6? 1 : atoi(argv[6]), atoi(argv[5]));

  return 0;
}
