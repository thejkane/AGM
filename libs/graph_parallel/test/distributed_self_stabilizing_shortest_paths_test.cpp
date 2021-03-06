// Copyright (C) 2005-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#define PBGL_ACCOUNTING

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/adjacency_list.hpp> // must precede delta_stepping...
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/distributed/self_stabilizing_shortest_paths.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <limits>
#include <iostream>

#include <am++/counter_coalesced_message_type.hpp>

#define CSR

#ifdef CSR
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#endif

/****************************************************************************
 * Edge weight generator iterator                                           *
 ****************************************************************************/
template<typename F, typename RandomGenerator>
class generator_iterator
{
public:
  typedef std::input_iterator_tag iterator_category;
  typedef typename F::result_type value_type;
  typedef const value_type&       reference;
  typedef const value_type*       pointer;
  typedef void                    difference_type;

  explicit 
  generator_iterator(RandomGenerator& gen, const F& f = F()) 
    : f(f), gen(&gen) 
  { 
    value = this->f(gen); 
  }

  reference operator*() const  { return value; }
  pointer   operator->() const { return &value; }

  generator_iterator& operator++()
  {
    value = f(*gen);
    return *this;
  }

  generator_iterator operator++(int)
  {
    generator_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const generator_iterator& other) const
  { return f == other.f; }

  bool operator!=(const generator_iterator& other) const
  { return !(*this == other); }

private:
  F f;
  RandomGenerator* gen;
  value_type value;
};

template<typename F, typename RandomGenerator>
inline generator_iterator<F, RandomGenerator> 
make_generator_iterator( RandomGenerator& gen, const F& f)
{ return generator_iterator<F, RandomGenerator>(gen, f); }

/****************************************************************************
 * Verification                                                             *
 ****************************************************************************/
template <typename Graph, typename DistanceMap, typename WeightMap>
void
verify_shortest_paths(Graph& g, DistanceMap distance, 
                      const WeightMap& weight) {

  distance.set_max_ghost_cells(0);

  {
    amplusplus::scoped_epoch epoch(g.transport());

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	get(distance, target(e, g));
      }
    }
  }

#ifdef PRINT_DEBUG
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
      if(get(distance, target(e, g)) > get(distance, source(e, g)) + get(weight, e))
	std::cerr << "ERROR: edge(" << source(e,g) << "," << target(e,g) << ") with distances (" << get(distance, source(e, g)) << "," << get(distance, target(e, g)) << "). Edge weight is " << get(weight, e) << ". Distance is " << get(distance, target(e, g)) - (get(distance, source(e, g)) + get(weight, e)) << " too long." << std::endl;
    }
  }
#endif
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
      assert(get(distance, target(e, g)) <= 
             get(distance, source(e, g)) + get(weight, e));
    }
  }
}

using namespace boost;

typedef int32_t weight_type; // Needed so atomics are available on MacOS

struct WeightedEdge {
  WeightedEdge(weight_type weight = 0) : weight(weight) { }
  
  weight_type weight;
};

struct VertexProperties {
  VertexProperties(int d = 0) : distance(d) { }

  weight_type distance;
};

namespace amplusplus {
  template<>
  struct make_mpi_datatype<WeightedEdge> : make_mpi_datatype_base {
    make_mpi_datatype<weight_type> dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1() {
      int blocklengths[1] = {1};
      MPI_Aint displacements[1];
      WeightedEdge test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.weight, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Datatype types[1] = {dt1.get()};
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };

  template<>
  struct make_mpi_datatype<VertexProperties> : make_mpi_datatype_base {
    make_mpi_datatype<weight_type> dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1() {
      int blocklengths[1] = {1};
      MPI_Aint displacements[1];
      VertexProperties test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.distance, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Datatype types[1] = {dt1.get()};
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}

template <typename Transport>
void 
test_distributed_shortest_paths(Transport& trans, int n, double p, int c, double test_p, double step_p, int seed, int num_threads)
{

#ifndef CSR
  typedef adjacency_list<vecS, 
                         distributedS<vecS>,
                         directedS,
                         VertexProperties,
                         WeightedEdge> Graph;
#else
  typedef compressed_sparse_row_graph<directedS, VertexProperties, WeightedEdge, no_property,
                                      distributedS<unsigned long> > Graph;
#endif

  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef property_map<Graph, vertex_index_t>::type vertex_index_map;
  typedef vector_property_map<vertex_descriptor, vertex_index_map> predecessor_map;

  // Build a random number generator
  minstd_rand gen;
  gen.seed(seed);

  // Build a random graph
#ifndef CSR
  Graph g(erdos_renyi_iterator<minstd_rand, Graph>(gen, n, p),
          erdos_renyi_iterator<minstd_rand, Graph>(),
          make_generator_iterator(gen, uniform_int<weight_type>(1, c)),
          n, trans);
#else
  Graph g(edges_are_unsorted, erdos_renyi_iterator<minstd_rand, Graph>(gen, n, p), erdos_renyi_iterator<minstd_rand, Graph>(), make_generator_iterator(gen, uniform_int<weight_type>(1, c)), n, trans);

#ifdef PRINT_DEBUG
  BGL_FORALL_VERTICES(v, g, Graph) { 
    BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
      std::cout << "(" << v << "," << target(e, g) << ")\n";
    }
  }
  std::cout << std::flush;
#endif

#endif



  // Choose a source vertex
  uniform_int<vertices_size_type> rand_vertex(0, n);
  vertex_descriptor source = vertex(rand_vertex(gen), g);

  // Make sure g is connected
  typedef two_bit_color_map<property_map<Graph, vertex_index_t>::type> ColorMap;
  ColorMap color(num_vertices(g), get(vertex_index, g));

  boost::graph::distributed::breadth_first_search<Graph> 
    bfs(g, make_bfs_visitor(null_visitor()), color);
  bfs.set_source(source);

  bfs(0);

  // predecessor_map predecessor(num_vertices(g), get(vertex_index, g));
  // vector_property_map<graph_traits<Graph>::vertex_descriptor, property_map<Graph, vertex_index_t>::type> predecessor(num_vertices(g), get(vertex_index, g));
  std::vector<vertex_descriptor> parentS(num_vertices(g), graph_traits<Graph>::null_vertex());
  iterator_property_map<std::vector<vertex_descriptor>::iterator, vertex_index_map> predecessor(parentS.begin(), get(vertex_index, g));

  typedef property_traits<ColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;
  BGL_FORALL_VERTICES(v, g, Graph) { 
    if (get(color, v) != Color::black()) {
      std::cout << "Graph must be connected\n";
      exit(-1);
    }
    //put(get(&VertexProperties::distance, g), v, std::numeric_limits<weight_type>::max());
    // put(predecessor, v, v); // Initialize to self
  }

  // If g is connected run delta-stepping and verify
  boost::graph::distributed::self_stabilizing_shortest_paths<Graph, 
      typename boost::property_map<Graph, weight_type VertexProperties::*>::type, 
      typename boost::property_map<Graph, weight_type WeightedEdge::*>::type,
							     iterator_property_map<std::vector<vertex_descriptor>::iterator, vertex_index_map> >
							     // predecessor_map>
    D(g, get(&VertexProperties::distance, g), get(&WeightedEdge::weight, g), predecessor, test_p, step_p);

  D.set_source(source);

  // Many threads now
  trans.set_nthreads(num_threads);

  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads]);
  for (int i = 0; i < num_threads; ++i) {
    std::cout << "starting thread" << i << std::endl;
    boost::thread thr(boost::ref(D), i);
    threads[i].swap(thr);
  }

  for (int i = 0; i < num_threads; ++i) {
    threads[i].join();
    std::cout << "finishing thread " << i << std::endl;
  }    

  // Back to one thread
  trans.set_nthreads(1);

  std::cout << "Finished all threads. Now verifying." << std::endl;

  verify_shortest_paths(g, 
                        get(&VertexProperties::distance, g), 
                        get(&WeightedEdge::weight, g));
}

int test_main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  int n = 1000;
  double p = 0.01;
  int c = 100;
  int seed = 1;
  int num_threads = 2;
  double test_p = 0.001;
  double step_p = 0.01;

  if (argc > 1) n             = lexical_cast<int>(argv[1]);
  if (argc > 2) p             = lexical_cast<double>(argv[2]);
  if (argc > 3) c             = lexical_cast<int>(argv[3]);
  if (argc > 4) test_p        = lexical_cast<double>(argv[4]);
  if (argc > 5) step_p        = lexical_cast<double>(argv[5]);
  if (argc > 6) seed          = lexical_cast<int>(argv[6]);
  if (argc > 7) num_threads   = lexical_cast<int>(argv[7]);
						
  test_distributed_shortest_paths(trans, n, p, c, test_p, step_p, seed, num_threads);

  return 0;
}
