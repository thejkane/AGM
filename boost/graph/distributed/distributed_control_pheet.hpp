// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_GRAPH_DISTRIBUTED_CONTROL_PHEET
#define BOOST_GRAPH_DISTRIBUTED_CONTROL_PHEET

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
#include <boost/graph/distributed/priority_q_defs.hpp>

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
	 //typename ProcessedInEdgeCountMap,
	 typename PriorityQueueGenerator = pheet_priority_queue_gen_128,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class distributed_control_pheet {
  // Type definitions for Vertex and Distance
  // Vertex information is in graph_traits::vertex_descriptor
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  // The distance is in the property traits of the edge
  typedef typename property_traits<EdgeWeightMap>::value_type Distance;
  // This is going to be our message
  typedef std::pair<Vertex, Distance> vertex_distance_data;

  // Default comparer
  struct default_comparer {
    bool operator()(const vertex_distance_data& vd1, const vertex_distance_data& vd2) {
      return vd1.second < vd2.second;
    }
  };

  // Owner map type - To calculate which vertex belong to which node
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  // Create the default queue generator type
  typedef typename PriorityQueueGenerator::template queue<vertex_distance_data, default_comparer>::type PriorityQueueType;

  // The handler, forward decleration
  struct vertex_distance_handler;

  typedef typename MessageGenerator::template call_result<
    vertex_distance_data,
    vertex_distance_handler,
    vertex_distance_owner<OwnerMap, vertex_distance_data>,
    amplusplus::idempotent_combination_t<boost::parallel::minimum<Distance>, Distance> >::type
  RelaxMessage;

public:

  distributed_control_pheet(Graph& g,
			   DistanceMap distance_map,
			   EdgeWeightMap weight_map,
			   amplusplus::transport &t,
			   MessageGenerator message_gen =
			   MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)),
			   MessageGenerator priority_message_gen =
			   MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12), 1),
			   unsigned int flushFreq = 2000,
			   unsigned int eagerLimit = 1000)
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_distance_data>(), 0)),
    graph(g),
    transport(t),
    distance_map(distance_map),
    weight_map(weight_map),
    //processed_in_edge_count_map(processed_in_edge_counts),
    owner_map(get(vertex_owner, g)),
    relax_msg(message_gen,
	      transport,
	      vertex_distance_owner<OwnerMap, vertex_distance_data>(owner_map),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Distance>(), std::numeric_limits<Distance>::max())),
    priority_relax_msg(priority_message_gen,
	      transport,
	      vertex_distance_owner<OwnerMap, vertex_distance_data>(owner_map),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Distance>(), std::numeric_limits<Distance>::max())),
    flushFrequency(flushFreq),
    eagerLimit(eagerLimit),
    pq(t.get_nthreads()),
    max_q_size(1000)
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

#ifdef PBGL2_PRINT_WORK_STATS
  // useful, invalidated, useless, rejected
  typedef std::tuple<unsigned long long, unsigned long long, unsigned long long, unsigned long long> work_stats_t;
  work_stats_t
  get_work_stats() {
    work_stats_t stats{ 0ul, 0ul, 0ul, 0ul };
    for(unsigned int i = 0; i < useful.size(); ++i) {
      std::get<0>(stats) += useful[i];
      std::get<1>(stats) += invalidated[i];
      std::get<2>(stats) += useless[i];
      std::get<3>(stats) += rejected[i];
    }
    return stats;
  }
#endif // PBGL2_PRINT_WORK_STATS

  size_t get_max_q_size() {
    return max_q_size;
  }

  void set_source(Vertex s) { source = s; }
  // What is actually getting called is this operator
  void operator() (int tid) { run(source, tid); }
  void run(Vertex s, int tid = 0);
  time_type get_start_time(){ return start_time; }
private:
  void initialize();
  void relax(const vertex_distance_data&, const int);
  void handle_queue(const int, amplusplus::transport::end_epoch_request&);

  const int dummy_first_member_for_init_order; // Unused
  const Graph& graph;
  amplusplus::transport& transport;
  DistanceMap distance_map;
  EdgeWeightMap weight_map;
  OwnerMap owner_map;
  Vertex source;
  unsigned int flushFrequency;
  unsigned int eagerLimit;
  RelaxMessage relax_msg;
  RelaxMessage priority_relax_msg;
  PriorityQueueType pq;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  unsigned int in_edge_factor;
  time_type start_time;
#ifdef DC_USE_PRIORITY
  std::vector<Distance> current_min;
#endif // DC_USE_PRIORITY
#ifdef PBGL2_PRINT_WORK_STATS
  std::vector<unsigned long long> useful, invalidated, useless, rejected;
#endif // PBGL2_PRINT_WORK_STATS
  size_t max_q_size;
};

#define DISTRIBUTE_CONTROL_PHEET_PARMS					\
  typename Graph, typename DistanceMap, typename EdgeWeightMap, typename PriorityQueueGenerator, typename MessageGenerator

#define DISTRIBUTE_CONTROL_PHEET_TYPE						\
  distributed_control_pheet<Graph, DistanceMap, EdgeWeightMap, PriorityQueueGenerator, MessageGenerator>

template<DISTRIBUTE_CONTROL_PHEET_PARMS> void
DISTRIBUTE_CONTROL_PHEET_TYPE::initialize() {
#ifdef PBGL2_PRINT_WORK_STATS
  useful.resize(transport.get_nthreads(), 0);
  invalidated.resize(transport.get_nthreads(), 0);
  useless.resize(transport.get_nthreads(), 0);
  rejected.resize(transport.get_nthreads(), 0);
#endif // PBGL2_PRINT_WORK_STATS
  relax_msg.set_handler(vertex_distance_handler(*this));
  priority_relax_msg.set_handler(vertex_distance_handler(*this));
}

template<DISTRIBUTE_CONTROL_PHEET_PARMS> void
DISTRIBUTE_CONTROL_PHEET_TYPE::run(Vertex s, int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    if (0 == tid) {
      int nthreads = transport.get_nthreads();

      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    }

    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized
    { amplusplus::scoped_epoch epoch(transport); }

    // Initialize priority queues
    pq.initialize(tid);

    t_bar->wait();

    if (tid == 0) {
      pq.re_group();
    }

    t_bar->wait();


    start_time = get_time();
    transport.begin_epoch();

    // Start the algorithm
    // If this is main thread and if node is the owner of the vertex do relax
    // for source with distance 0
    if (get(owner_map, s) == transport.rank() && tid == 0) {
#ifdef PRINT_DEBUG
      std::cout << "Increasing active count at source relax " << std::endl;
#endif
      transport.increase_activity_count(1);
      relax(vertex_distance_data(s, 0), tid);
    }

    // Check whether we still have work to do

    amplusplus::transport::end_epoch_request request = transport.i_end_epoch();
    // if we have work process queues. queue is populated in the relax_msg
    handle_queue(tid, request);
  }

}

template<DISTRIBUTE_CONTROL_PHEET_PARMS>
void DISTRIBUTE_CONTROL_PHEET_TYPE::handle_queue(const int tid, amplusplus::transport::end_epoch_request& request) {
  //int doFlushCounter = 0;

  std::pair<Vertex, Distance> vd;
  while (true) {
    if (pq.try_pop(vd)) {

      // get the top vertex
      Vertex v = vd.first;
      Distance v_distance = vd.second;

      // remember if the queue was empty after pop
      //      bool was_empty = pq.is_empty(tid);

      // get existing distance
      // if distance in distance map is better that distance we got from queue
      // then we dont need to do anything (i.e. another thread has updated distance map with better
      // distance), otherwise send new distances for neighbours
      Distance map_dist = get(distance_map, v);

      if (map_dist == v_distance) {
	BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
	  Vertex u = target(e, graph);
	  Distance we = get(weight_map, e);
	  // Increase active count
#ifdef PRINT_DEBUG
	  std::cout << "Increasing active count at message sending " << std::endl;
#endif
	  transport.increase_activity_count(1);
	  relax_msg.send(vertex_distance_data(u, v_distance + we));
	}
      }
#ifdef PBGL2_PRINT_WORK_STATS
      else {
        ++rejected[tid];
      }
#endif

#ifdef PRINT_DEBUG
      std::cout << "Decreasing activity count at the end of pop" << std::endl;
#endif
      //      completed_count++;
      transport.decrease_activity_count(1);
    } else {
      if (request.test()) break;
    }
  }
}

template<DISTRIBUTE_CONTROL_PHEET_PARMS> void
DISTRIBUTE_CONTROL_PHEET_TYPE::relax(const vertex_distance_data& vd, const int tid) {
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
      pq.push(vd);

#ifdef PBGL2_PRINT_WORK_STATS
      ++useful[tid];
      if(old_dist < std::numeric_limits<Distance>::max()) ++invalidated[tid];
#endif // PBGL2_PRINT_WORK_STATS
      return;
    } // end of if (last_old_dist == old_dist) {
  } // end of while (new_distance < old_dist) {
#ifdef PBGL2_PRINT_WORK_STATS
  ++useless[tid];
#endif // PBGL2_PRINT_WORK_STATS

#ifdef PRINT_DEBUG
  std::cout << "Decreasing activity count at end of relax" << std::endl;
#endif
  transport.decrease_activity_count(1);
}

template<DISTRIBUTE_CONTROL_PHEET_PARMS>
struct DISTRIBUTE_CONTROL_PHEET_TYPE::vertex_distance_handler {
  vertex_distance_handler() : self(NULL) {}
  vertex_distance_handler(distributed_control_pheet& self) : self(&self) {}

  void operator() (const vertex_distance_data& data) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->relax(data, tid);
  }

protected:
  distributed_control_pheet* self;
};

} // end namespace distributed
} // end namespace graph
} // end namespace boost

#endif
