// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_GRAPH_DISTRIBUTED_CONTROL
#define BOOST_GRAPH_DISTRIBUTED_CONTROL

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

// The default priority queue
template<typename V, typename D, typename Compare>
struct default_priority_queue {

  typedef std::priority_queue<std::pair<V, D>, std::vector<std::pair<V, D> >, Compare> DefaultPriorityQueueType;

  DefaultPriorityQueueType priority_q;

  void put(const std::pair<V, D>& p) {
    priority_q.push(p);
  }

  std::pair<V, D> pop() {
    // get the top vertex
    std::pair<V, D> topv = priority_q.top();
    // remove top element from the queue
    priority_q.pop();
    return topv;
  }

  bool is_empty() const {
    return priority_q.empty();
  }

  std::pair<V, D> top() {
    // get the top vertex
    std::pair<V, D> topv = priority_q.top();
    return topv;
  }

  size_t size() const {
    return priority_q.size();
  }

};

// Default priority queue generator
struct default_priority_queue_gen {
  template<typename V, typename D, typename Compare>
  struct queue {
    typedef default_priority_queue<V, D, Compare> type;
  };
};

// Vector of vector implementation
/*template<typename V, typename D, typename Compare>
class vector_of_vector_table {
public:
  vector_of_vector_table()
    : max_distance(0),min_distance(std::numeric_limits<int32_t>::max()),_size(0){}

  void put(const std::pair<V, D>& p) {
    //if the pair passed to has greater distance than the max, resize more
    //then push
    D current_distance = p.second;
#ifdef PRINT_DEBUG
    std::cerr<<"pushing vertex "<<p.first<< " distance processing "<<current_distance<<"\n";
#endif

    if (current_distance >= max_distance) {
      buckets.resize(current_distance+1);
      max_distance = current_distance;
    }
    if (current_distance <= min_distance) {
      min_distance = current_distance;
    }
    ++_size;
    if(buckets.size() < current_distance || buckets.size() == 0) {
      if(buckets.size() == 0){
	buckets.resize(1);
	//min_distance = 0;

      }
      else
	buckets.resize(current_distance+1);
    }
    buckets[current_distance].push_back(p.first);

#ifdef PRINT_DEBUG
    std::cerr<<"PUSH min: "<< min_distance<< " Max "<<max_distance<<"\n";
#endif

  }

  std::pair<V, D> pop() {
    // get the top vertex

    // remove top element from the bucket
    --_size;
    V& topv = buckets[min_distance].back();
    D distance_ret = min_distance;

#ifdef PRINT_DEBUG
    std::cerr<<"POP min"<<min_distance<<" \n";
    std::cerr<<"Popping "<< topv <<" with distance " << min_distance <<"\n";
#endif

    buckets[min_distance].pop_back();
    if(buckets[min_distance].empty()) {
      for(auto i = min_distance+1; i <= max_distance;i++) {
	if(!buckets[i].empty()){
	  min_distance = i;
	  break;
	}
      }
    }

    if (_size==0) {
      min_distance = std::numeric_limits<int32_t>::max();
      max_distance = 0;
    }
    return std::make_pair(topv, distance_ret);
  }

  bool is_empty() const {
    return (_size==0);
  }


  std::pair<V, D> top() {
    // get the top vertex
    V& topv = buckets[min_distance].back();
    D distance_ret = min_distance;
    return std::make_pair(topv, distance_ret);
  }
  size_t size() const {
    return _size;

  }

private:
  typedef std::vector<std::vector<V> > Buckets;

  Buckets buckets;
  D max_distance;
  D min_distance;
  int _size;
};


// vector of vector generator
struct vector_of_vector_gen {
  template<typename V, typename D, typename Compare>
  struct queue {
	typedef vector_of_vector_table<V, D, Compare> type;
  };

  };*/

template<typename Graph,
	 typename DistanceMap,
	 typename EdgeWeightMap,
	 typename WorkStats,
	 typename PriorityQueueGenerator = default_priority_queue_gen,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class distributed_control {
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

  // Create the default queue generator type
  typedef typename PriorityQueueGenerator::template queue<Vertex, Distance, default_comparer>::type PriorityQueueType;

  // The handler, forward decleration
  struct vertex_distance_handler;

  typedef typename MessageGenerator::template call_result<
    vertex_distance_data,
    vertex_distance_handler,
    vertex_distance_owner<OwnerMap, vertex_distance_data>,
    amplusplus::idempotent_combination_t<boost::parallel::minimum<Distance>, Distance> >::type
  RelaxMessage;

public:

  distributed_control(Graph& g,
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
  void relax(const vertex_distance_data&, const int);
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
  std::vector<PriorityQueueType> thread_buckets;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  unsigned int in_edge_factor;
  time_type start_time;
#ifdef DC_USE_PRIORITY
  std::vector<Distance> current_min;
#endif // DC_USE_PRIORITY
};

#define DISTRIBUTE_CONTROL_PARMS					\
  typename Graph, typename DistanceMap, typename EdgeWeightMap, typename WorkStats, typename PriorityQueueGenerator, typename MessageGenerator

#define DISTRIBUTE_CONTROL_TYPE						\
  distributed_control<Graph, DistanceMap, EdgeWeightMap, WorkStats, PriorityQueueGenerator, MessageGenerator>

template<DISTRIBUTE_CONTROL_PARMS> void
DISTRIBUTE_CONTROL_TYPE::initialize() {
  relax_msg.set_handler(vertex_distance_handler(*this));

  // Initialize distance labels
  BGL_FORALL_VERTICES_T(v, graph, Graph) { 
    put(distance_map, v, (std::numeric_limits<Distance>::max)()); 
  }
}

template<DISTRIBUTE_CONTROL_PARMS> void
DISTRIBUTE_CONTROL_TYPE::run(Vertex s, int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    if (0 == tid) {
      int nthreads = transport.get_nthreads();

      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));

      // Resize thread queue
      thread_buckets.resize(nthreads);

#ifdef DC_USE_PRIORITY
      // Resize current minimums
      current_min.resize(nthreads);
#endif // DC_USE_PRIORITY

      for (int k = 0; k < nthreads; ++k) {
	// For each thread create a priority queue
	thread_buckets[k] = PriorityQueueType();
#ifdef DC_USE_PRIORITY
	// For each thread set the minimum to 0
	current_min[k] = 0;
#endif // DC_USE_PRIORITY
      }
    }

    // TODO delta had following. Not sure we also need this - ask
    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized
    { amplusplus::scoped_epoch epoch(transport); }

    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
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

    // should come before begin epoch
    start_time = get_time();
    // Start the algorithm
    transport.begin_epoch();
    // If this is main thread and if node is the owner of the vertex do relax
    // for source with distance 0
    if (get(owner_map, s) == transport.rank() && tid == 0) {
      relax(vertex_distance_data(s, 0), tid);
    }

    // Check whether we still have work to do

    amplusplus::transport::end_epoch_request request = transport.i_end_epoch();
    // if we have work process queues. queue is populated in the relax_msg
    handle_queue(tid, request);
  }

#ifdef PRINT_DEBUG
  std::cerr<<"Done with run\n";
#endif
}

template<DISTRIBUTE_CONTROL_PARMS>
void DISTRIBUTE_CONTROL_TYPE::handle_queue(const int tid, amplusplus::transport::end_epoch_request& request) {
  // get the queue
  PriorityQueueType& queue = thread_buckets[tid];
  int doFlushCounter = 0;


  while (true) {
    if (!queue.is_empty()) {
      // get the top vertex
      std::pair<Vertex, Distance> vd = queue.pop();
      Vertex v = vd.first;
      Distance v_distance = vd.second;

#ifdef DC_USE_PRIORITY
      current_min[tid] = vd.second;
#endif // DC_USE_PRIORITY

      // remember if the queue was empty after pop
      bool was_empty = queue.is_empty();

      // get existing distance
      // if distance in distance map is better that distance we got from queue
      // then we dont need to do anything (i.e. another thread has updated distance map with better
      // distance), otherwise send new distances for neighbours
      Distance map_dist = get(distance_map, v);

      if (map_dist == v_distance) {
	BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
	  Vertex u = target(e, graph);
	  Distance we = get(weight_map, e);
	  relax_msg.send(vertex_distance_data(u, v_distance + we));
	}
      }

      // If the queue was empty after the pop, we need to decrease activity count.  The queue may not be empty at this point, because a send could result in a push to the queue.
      if (was_empty) {

#ifdef PRINT_DEBUG
	std::cout << "########### Decrease activity count by 1 ########## " << std::endl;
#endif
	// after poping element if queue is empty decrease the activity count
	transport.decrease_activity_count(1);
      }

      doFlushCounter++;

      if(doFlushCounter == flushFrequency) {
	doFlushCounter = 0;

	if(request.test()) {
#ifdef PRINT_DEBUG
	  assert(queue.is_empty());
#endif
	  break;
	}
      }

    } else {
      if (request.test()) break;
    }
  }
}

template<DISTRIBUTE_CONTROL_PARMS> void
DISTRIBUTE_CONTROL_TYPE::relax(const vertex_distance_data& vd, const int tid) {

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
      PriorityQueueType& queue = thread_buckets[tid];

      if (queue.is_empty()) {
	// If the queue was empty, increase activity count to notify AM++ that work is being held off.
#ifdef PRINT_DEBUG
        std::cout << "########### Increasing activity count by 1 ########## " << std::endl;
#endif
        transport.increase_activity_count(1);
      }

      queue.put(vd);

#ifdef PBGL2_PRINT_WORK_STATS
      int tid = amplusplus::detail::get_thread_id();
      if(old_dist < std::numeric_limits<Distance>::max()) { 
	work_stats.increment_invalidated(tid);
      } else {
	work_stats.increment_useful(tid);
      }
#endif // PBGL2_PRINT_WORK_STATS

      return;
    }
  }
#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(amplusplus::detail::get_thread_id());
#endif // PBGL2_PRINT_WORK_STATS
}

template<DISTRIBUTE_CONTROL_PARMS>
struct DISTRIBUTE_CONTROL_TYPE::vertex_distance_handler {
  vertex_distance_handler() : self(NULL) {}
  vertex_distance_handler(distributed_control& self) : self(&self) {}

  void operator() (const vertex_distance_data& data) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->relax(data, tid);
  }

protected:
  distributed_control* self;
};

} // end namespace distributed
} // end namespace graph
} // end namespace boost

#endif
