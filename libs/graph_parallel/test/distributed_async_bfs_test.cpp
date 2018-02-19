// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds

// Simple asyncrhonous BFS implementation, of course you could also
// use async_breadth_first_search from
// boost/graph/distributed/breadth_first_searc.hpp

#include <am++/am++.hpp> 
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/test/minimal.hpp>
#include <boost/lexical_cast.hpp>

#define COALESCED_SIZE 4096

using namespace boost;

typedef boost::graph_traits<compressed_sparse_row_graph<directedS, 
			    no_property, no_property, no_property, 
			    distributedS<> > >::vertex_descriptor Vertex;
typedef std::pair<Vertex, std::size_t> bfs_distance_data;

template <typename Graph, typename DistanceMap>
struct async_bfs_handler {

  typedef amplusplus::basic_coalesced_message_type<bfs_distance_data, async_bfs_handler>
  async_bfs_message_type;

  typedef typename property_traits<DistanceMap>::value_type distance_type;

  async_bfs_handler() {}
  async_bfs_handler(Graph& g, DistanceMap& dist, async_bfs_message_type& msg) 
    : g(&g), distance(&dist), bfs_msg(&msg) {}
  
  void operator() (int /* source */, const bfs_distance_data& data)
  {
    typename graph_traits<Graph>::vertex_descriptor v = data.first;
    /* volatile */ distance_type old_dist;
    distance_type new_dist = data.second;

    do {
      old_dist = get(*distance, v);
      if (new_dist >= old_dist) return;
    } while (!exchange(*distance, v, old_dist, new_dist));

    // TODO: Explore all out edges of v
    BGL_FORALL_ADJ_T(v, u, *g, Graph) { 
      bfs_msg->send(std::make_pair(u, new_dist + 1), get(get(vertex_owner, *g), u)); 
    }
  }
  
protected:
  Graph*                  g;
  DistanceMap*            distance;
  async_bfs_message_type* bfs_msg;
};

int test_main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  int num_threads = 1;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "--threads")
      num_threads = boost::lexical_cast<int>( argv[i+1] );
  }

  typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property, 
                                      distributedS<> >
    Digraph;

  // Build a Watts-Strogatz graph to test with
  typedef small_world_iterator<minstd_rand, Digraph> SWIter;

  int n = 100, k = 4;
  double prob = 0.05;

  minstd_rand gen;
  gen.seed(1);
  Digraph g(edges_are_unsorted, SWIter(gen, n, k, prob), SWIter(), n, trans);

  // Distance map
  std::vector<std::size_t> distanceS(num_vertices(g), std::numeric_limits<std::size_t>::max());
  typedef iterator_property_map<std::vector<std::size_t>::iterator, 
                                property_map<Digraph, vertex_index_t>::type> 
    DistanceMap;
  DistanceMap distance(distanceS.begin(), get(vertex_index, g));

  // Messages
  amplusplus::register_mpi_datatype<bfs_distance_data>();
  typedef async_bfs_handler<Digraph, DistanceMap> Handler;
  typedef amplusplus::basic_coalesced_message_type<bfs_distance_data, Handler>
    async_bfs_message_type;
  async_bfs_message_type async_bfs_msg(amplusplus::basic_coalesced_message_type_gen(COALESCED_SIZE), trans);

  //
  // TODO: Add addressing to message type
  //

  async_bfs_msg.set_handler(Handler(g, distance, async_bfs_msg));

  // Single threaded since we'd have to fork inside the coalesced messages
  {
    amplusplus::scoped_epoch epoch(trans);
    async_bfs_msg.send(std::make_pair(vertex(0, g), 0), get(get(vertex_owner, g), vertex(0, g)));
  }

  int visited = 0;
  BGL_FORALL_VERTICES(v, g, Digraph) { 
    if (get(distance, v) < std::numeric_limits<std::size_t>::max())
      ++visited; 
  }

  std::cerr << trans.rank() << ": visited " << visited << " vertices of " << num_vertices(g) << std::endl;

  return 0;
}
