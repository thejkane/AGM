// Copyright (C) 2011 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds

// A test of the distributed breadth first search algorithm

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/test/minimal.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/parallel/algorithm.hpp>

#include <queue>

using namespace boost;

/****************************************************************************
 * Edge Property                                                            *
 ****************************************************************************/
typedef int weight_type;

struct WeightedEdge {
  WeightedEdge(weight_type weight = 0) : weight(weight) { }
  
  weight_type weight;
};

namespace amplusplus {
  template<>
  struct make_mpi_datatype<WeightedEdge> : make_mpi_datatype_base {
    make_mpi_datatype<weight_type> dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype() : dt1() {
      int blocklengths[1] = {1};
      MPI_Aint displacements[1];
      char dummy;
      WeightedEdge *test_object = (WeightedEdge*)(&dummy);
      MPI_Aint test_object_ptr;
      MPI_Get_address(test_object, &test_object_ptr);
      MPI_Get_address(&test_object->weight, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Datatype types[1] = { dt1.get() };
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}

typedef graph_traits<compressed_sparse_row_graph<directedS, no_property, 
  WeightedEdge, no_property, distributedS<> > >::vertex_descriptor 
    Vertex;

typedef boost::tuple<Vertex, weight_type> vertex_distance_data;

template <typename Compare>
struct compare_second {

  inline compare_second(const Compare& c = Compare()) : cmp(c) { }

  template <typename A, typename B>
  inline bool operator()(const boost::tuple<A, B>& x, const boost::tuple<A, B>& y) const
  {
    return cmp(get<1>(x), get<1>(y));
  }

  protected:
    Compare cmp;
};

template <typename Graph, typename PMap>
struct get_pmap_of_first {

  typedef typename property_traits<PMap>::value_type value_type;
  typedef typename property_traits<PMap>::reference reference;
  typedef typename property_traits<PMap>::category category;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef vertex_distance_data key_type;

  inline get_pmap_of_first(const PMap& p) : p(p) {}

  value_type get(const key_type& key) 
  { return boost::get(p, tuples::get<0>(key)); }

private:
  const PMap& p;
};

template <typename Graph, typename PMap>
typename property_traits<PMap>::value_type
get(get_pmap_of_first<Graph, PMap> pmap, 
    const typename get_pmap_of_first<Graph, PMap>::key_type& key)
{ return pmap.get(key); }

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

template <typename T, typename Sequence, typename Compare>
class my_priority_queue 
  : public std::priority_queue<T, Sequence, Compare>
{
  typedef typename std::priority_queue<T, Sequence, Compare> inherited;

public:
  void push_back(const typename inherited::value_type& v) { inherited::push(v); }
};

template <typename Graph, typename WeightMap, typename DistanceMap, 
	  typename Queue>
void 
run(Graph& g, WeightMap weight, DistanceMap distance, Queue& Q,
    typename graph_traits<Graph>::vertex_descriptor source, 
    weight_type lookahead)
{
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;

  // Push source
  {
    amplusplus::scoped_epoch epoch(g.transport());

    if (get(get(vertex_owner, g), source) == g.transport().rank())
      Q.push(boost::make_tuple(source, weight_type(0)));
  }

  // Run dijkstra
  weight_type min_distance = 0;
  Vertex v;
  weight_type dist;
  std::size_t queue_size, global_queue_size;
  std::size_t epoch_count = 0;

  do {
    
    ++epoch_count;
    if (epoch_count > 10) 
      lookahead *= 10; 

    std::cout << g.transport().rank() << ": EPOCH START\n";

    std::cout << g.transport().rank() << ": Local queue size: " << Q.size() << std::endl; 
    
    queue_size = Q.size();

    weight_type local_min = Q.empty() ? std::numeric_limits<weight_type>::max() : get<1>(Q.top());

    using boost::parallel::all_reduce;
    using boost::parallel::minimum;

    all_reduce<weight_type, minimum<weight_type> > r(g.transport(), minimum<weight_type>());
    min_distance = r(local_min);

    if (g.transport().rank() == 0)
      std::cout << "Global min distance: " << min_distance << std::endl;

    {
      amplusplus::scoped_epoch_value(g.transport(), queue_size, global_queue_size);

      while(!Q.empty() && get<1>(Q.top()) <= min_distance + lookahead) {
	tuples::tie(v, dist) = Q.top(); Q.pop();
	if (dist < get(distance, v)) {
	  put(distance, v, dist);
	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    // We could filter local edges that don't need to be queued here
	    Q.push(boost::make_tuple(target(e, g), get(distance, v) + get(weight, e))); 
	  }
	}
      }
    }

    if (g.transport().rank() == 0)
      std::cout << "Global queue size: " << global_queue_size << std::endl;

    std::cout << g.transport().rank() << ": EPOCH END\n";    

  } while(global_queue_size > 0);
}

int test_main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge, 
                                      no_property, distributedS<> >
      Digraph;
  
  typedef graph_traits<Digraph>::vertex_descriptor Vertex;
  typedef property_map<Digraph, vertex_index_t>::type VertexIndexMap;
  typedef property_map<Digraph, vertex_owner_t>::type OwnerMap;
  
  // Build an Erdos-Renyi graph to test with
#ifdef ER
  typedef erdos_renyi_iterator<minstd_rand, Digraph> ERIter;
#else
  typedef rmat_iterator<minstd_rand, Digraph> RMATIter;
#endif  

  int n = 100, C = 10, m = 500;
  double prob;
  int source_idx = 0;
  weight_type lookahead = 1.5*C;

  // Parse args
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "--vertices")
      n = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--max-edge-weight")
      C = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--edges")
      m = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--prob")
      prob = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--source")
      source_idx = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--lookahead")
      lookahead = boost::lexical_cast<size_t>( argv[i+1] );
  }

  prob = static_cast<double>(m) / n / n;

  minstd_rand gen;
  gen.seed(1);

#ifdef ER
  Digraph g(edges_are_unsorted, ERIter(gen, n, prob), ERIter(),
#else
  Digraph g(edges_are_unsorted, RMATIter(gen, n, m, 0.57, 0.19, 0.19, 0.05), RMATIter(), 
#endif
	    make_generator_iterator(gen, uniform_int<int>(1, C)), 
	    n, trans);

  // Output edge list  
#ifdef OUTPUT_GRAPH
  BGL_FORALL_EDGES(e, g, Digraph) {
    std::cout << get(get(vertex_local, g), source(e, g)) << "@" 
	      << get(get(vertex_owner, g), source(e, g))
	      << "->" << get(get(vertex_local, g), target(e, g)) 
	      << "@" << get(get(vertex_owner, g), target(e, g)) << std::endl;
  }
#endif

  // Distance map
  std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
  typedef iterator_property_map<std::vector<weight_type>::iterator, 
                                property_map<Digraph, vertex_index_t>::type>
    DistanceMap;
  DistanceMap distance(distanceS.begin(), get(vertex_index, g));

  // Queue
  typedef vertex_distance_data queue_value_type;
  typedef compare_second<std::less<weight_type> > queue_compare_type;
  typedef get_pmap_of_first<Digraph, OwnerMap> IndirectOwnerMap;

  typedef my_priority_queue<queue_value_type, std::vector<queue_value_type>, queue_compare_type>
    local_queue_type;

  local_queue_type local_queue;

  amplusplus::register_mpi_datatype<queue_value_type>();
  boost::graph::distributed::distributed_queue<IndirectOwnerMap, local_queue_type>
	    Q(g.transport(), IndirectOwnerMap(get(vertex_owner, g)), boost::make_shared<local_queue_type>(local_queue));

  run(g, get(&WeightedEdge::weight, g), distance, Q, vertex(source_idx, g), lookahead);

  return 0;	    
}
