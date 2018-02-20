// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_DISTRIBUTED_CONTROL_CHAOTIC
#define BOOST_GRAPH_DISTRIBUTED_CONTROL_CHAOTIC

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/transport.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap

#include <queue>
#include <algorithm> // for std::min, std::max
#include <iostream>
#include <atomic>
#include <tuple>
#include <climits>

namespace boost {
namespace graph {
namespace distributed {

template<typename Graph,
	 typename DistanceMap,
	 typename EdgeWeightMap,
	 typename WorkStats,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class distributed_control_chaotic {
  // Type definitions for Vertex and Distance
  // Vertex information is in graph_traits::vertex_descriptor
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  // The distance is in the property traits of the edge
  typedef typename property_traits<EdgeWeightMap>::value_type Distance;
  // This is going to be our message
  typedef std::pair<Vertex, Distance> vertex_distance_data;

  // Default comparer
  // TODO: MZ: We can probably turn of sorting as discussed.  Sorting will occur in the priority queue.  Once we want to try intra-buffer sorting we can revisit this. - this can be done using the preprocessor macro. I guess once we integrate Jesun's code on command line switch we dont need the preprocessor macro
  struct default_comparer {
    bool operator()(const vertex_distance_data& vd1, const vertex_distance_data& vd2) {
      return vd1.second > vd2.second;
    }
  };

  // Owner map type - To calculate which vertex belong to which node
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  // The handler, forward decleration
  struct vertex_distance_handler;

  typedef typename MessageGenerator::template call_result<
    vertex_distance_data,
    vertex_distance_handler,
    vertex_distance_owner<OwnerMap, vertex_distance_data>,
    amplusplus::idempotent_combination_t<boost::parallel::minimum<Distance>, Distance> >::type
  RelaxMessage;


public:

  distributed_control_chaotic(Graph& g,
                              DistanceMap distance_map,
                              EdgeWeightMap weight_map,
                              amplusplus::transport &t,
                              int offs,
                              WorkStats& stats,
                              MessageGenerator message_gen =
                              MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)),
                              unsigned int flushFreq = 2000)
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_distance_data>(), 0)),
    graph(g),
    transport(t),
    core_offset(offs),
    work_stats(stats),
    distance_map(distance_map),
    weight_map(weight_map),
    owner_map(get(vertex_owner, g)),
    relax_msg(message_gen,
	      transport,
	      vertex_distance_owner<OwnerMap, vertex_distance_data>(owner_map),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Distance>(), std::numeric_limits<Distance>::max())),
    flushFrequency(flushFreq)
  {
    initialize();
  }


#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  typedef std::pair<unsigned long long, unsigned long long> cache_stats;
  cache_stats
  get_cache_stats() {
    cache_stats stats(0, 0);
    for(size_t i = 0; i < relax_msg.counters.hits.size(); ++i) {
      stats.first += relax_msg.counters.hits[i];
      stats.second += relax_msg.counters.tests[i];
    }
    return stats;
  }
#endif // AMPLUSPLUS_PRINT_HIT_RATES

  void set_source(Vertex s) { source = s; }
  // What is actually getting called is this operator
  void operator() (int tid) { run(source, tid); }
  void run(Vertex s, int tid = 0);
  time_type get_start_time() { return start_time; }
  /*void print_in_edge_map() {
    typename graph_traits<Graph>::vertex_iterator i, end;
    unsigned int k = 0;
    //    std::pair<vertex_iter, vertex_iter> vp;
    //for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
    for(tie(i, end) = vertices(graph); i != end; ++i) {
      ++k;
      unsigned int val = processed_in_edge_count_map[*i];
      std::cout << "The vertex " << k << " is visited " << val << " times." << std::endl;
      if (k == 100000)
	break;
    }
    }*/
private:
  void initialize();
  void relax(const vertex_distance_data&);
  void handle_queue(const int, amplusplus::transport::end_epoch_request&);

  const int dummy_first_member_for_init_order; // Unused
  const Graph& graph;
  amplusplus::transport& transport;
  int core_offset;
  WorkStats& work_stats;
  DistanceMap distance_map;
  EdgeWeightMap weight_map;
  //ProcessedInEdgeCountMap processed_in_edge_count_map;
  OwnerMap owner_map;
  Vertex source;
  unsigned int flushFrequency;
  RelaxMessage relax_msg;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  unsigned int in_edge_factor;
  time_type start_time;
};

#define DISTRIBUTE_CONTROL_CHAOTIC_PARMS					\
  typename Graph, typename DistanceMap, typename EdgeWeightMap, typename WorkStats, typename MessageGenerator

#define DISTRIBUTE_CONTROL_CHAOTIC_TYPE						\
  distributed_control_chaotic<Graph, DistanceMap, EdgeWeightMap, WorkStats, MessageGenerator>

template<DISTRIBUTE_CONTROL_CHAOTIC_PARMS> void
DISTRIBUTE_CONTROL_CHAOTIC_TYPE::initialize() {
  relax_msg.set_handler(vertex_distance_handler(*this));

  // Initialize distance labels
  BGL_FORALL_VERTICES_T(v, graph, Graph) { 
    put(distance_map, v, (std::numeric_limits<Distance>::max)()); 
  }
}

template<DISTRIBUTE_CONTROL_CHAOTIC_PARMS> void
DISTRIBUTE_CONTROL_CHAOTIC_TYPE::run(Vertex s, int tid) {

  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    if (0 == tid) {
      int nthreads = transport.get_nthreads();

      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    }
  
    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized
    { amplusplus::scoped_epoch epoch(transport); }

    t_bar->wait();
    if (pin(tid+core_offset) != 0) {
      std::cerr << "[ERROR] Unable to pin current thread to "
		<< "core : " << tid << std::endl;
      assert(false);
    }

    // wait till all threads are pinned
    t_bar->wait();
    { amplusplus::scoped_epoch epoch(transport); }

    validate_thread_core_relation();

    start_time = get_time();
    // Start the algorithm
    {
      amplusplus::scoped_epoch epoch(transport); 
      // If this is main thread and if node is the owner of the vertex do relax
      // for source with distance 0
      if (get(owner_map, s) == transport.rank() && tid == 0) {
	relax(vertex_distance_data(s, 0));
      }
    }

  }
#ifdef PRINT_DEBUG
  std::cerr<<"Done with run\n";
#endif
}

template<DISTRIBUTE_CONTROL_CHAOTIC_PARMS> void
DISTRIBUTE_CONTROL_CHAOTIC_TYPE::relax(const vertex_distance_data& vd) {

#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_edges(amplusplus::detail::get_thread_id());
#endif

  using boost::parallel::val_compare_and_swap;

  Vertex v = vd.first;
  Distance new_distance = vd.second;

  // we are processing an incoming edge for vertex v
  // increase incoming edge count
  //++processed_in_edge_count_map[v];

  // get existing distance
  // need last_old_dist to atomically update distance map
  Distance old_dist = get(distance_map, v), last_old_dist;

  // check whether we got a better distance if we got better
  // distance update the distance map
  // credit : took from delta_stepping_shortest_paths.hpp
  while (new_distance < old_dist) {
    last_old_dist = old_dist;
    old_dist = val_compare_and_swap(&distance_map[v], old_dist, new_distance);
    if (last_old_dist == old_dist) {
      // we got a better distance - insert it to work queue
#ifdef PBGL2_PRINT_WORK_STATS
      int tid = amplusplus::detail::get_thread_id();
      if(old_dist < std::numeric_limits<Distance>::max()) { 
	work_stats.increment_invalidated(tid);
      } else {
	work_stats.increment_useful(tid);
      }
#endif
      BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
	Vertex u = target(e, graph);
	Distance we = get(weight_map, e);
	relax_msg.send(vertex_distance_data(u, new_distance + we));
      }

      break;
    }
  }
  
#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(amplusplus::detail::get_thread_id());
#endif
}

template<DISTRIBUTE_CONTROL_CHAOTIC_PARMS>
struct DISTRIBUTE_CONTROL_CHAOTIC_TYPE::vertex_distance_handler {
  vertex_distance_handler() : self(NULL) {}
  vertex_distance_handler(distributed_control_chaotic& self) : self(&self) {}

  void operator() (const vertex_distance_data& data) const {
    self->relax(data);
  }

protected:
  distributed_control_chaotic* self;
};

} // end namespace distributed
} // end namespace graph
} // end namespace boost

#endif