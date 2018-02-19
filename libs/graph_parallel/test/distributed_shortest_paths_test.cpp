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
#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <boost/parallel/append_buffer.hpp>

#include <am++/counter_coalesced_message_type.hpp>

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

  (void)weight; // Avoid warning for unused weight with NDEBUG

  distance.set_max_ghost_cells(0);

  {
    amplusplus::scoped_epoch epoch(g.transport());

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	get(distance, target(e, g));
      }
    }
  }

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
test_distributed_shortest_paths(Transport& trans, int n, double p, int c, int seed, int num_threads)
{
  typedef adjacency_list<vecS, 
                         distributedS<>,
                         directedS,
                         VertexProperties,
                         WeightedEdge> Graph;

  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef property_map<Graph, vertex_index_t>::type vertex_index_map;

  // Build a random number generator
  minstd_rand gen;
  gen.seed(seed);

  // Build a random graph
  Graph g(erdos_renyi_iterator<minstd_rand, Graph>(gen, n, p),
          erdos_renyi_iterator<minstd_rand, Graph>(),
          make_generator_iterator(gen, uniform_int<weight_type>(1, c)),
          n, trans);

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

  typedef property_traits<ColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;
  BGL_FORALL_VERTICES(v, g, Graph) { 
    if (get(color, v) != Color::black()) {
      std::cout << "Graph must be connected\n";
      exit(-1);
    }
  }

  // If g is connected run delta-stepping and verify
  boost::graph::distributed::delta_stepping_shortest_paths<Graph, 
      typename boost::property_map<Graph, weight_type VertexProperties::*>::type, 
      typename boost::property_map<Graph, weight_type WeightedEdge::*>::type> 
    D(g, get(&VertexProperties::distance, g), get(&WeightedEdge::weight, g));

  D.set_source(source);

  // Many threads now
  trans.set_nthreads(num_threads);

  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads]);
  for (int i = 0; i < num_threads; ++i) {
    boost::thread thr(boost::ref(D), i);
    threads[i].swap(thr);
  }

  for (int i = 0; i < num_threads; ++i)
    threads[i].join();
    
  // Back to one thread
  trans.set_nthreads(1);

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

  if (argc > 1) n             = lexical_cast<int>(argv[1]);
  if (argc > 2) p             = lexical_cast<double>(argv[2]);
  if (argc > 3) c             = lexical_cast<int>(argv[3]);
  if (argc > 4) seed          = lexical_cast<int>(argv[4]);
  if (argc > 5) num_threads   = lexical_cast<int>(argv[5]);

  test_distributed_shortest_paths(trans, n, p, c, seed, num_threads);

  return 0;
}
