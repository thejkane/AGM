// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Connected Components Algortihm================//
// Distributed control version of connected components algorithm.
//===========================================================//

#ifndef DATA_DRIVEN_CONNECTED_COMPONENTS
#define DATA_DRIVEN_CONNECTED_COMPONENTS
#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/transport.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <boost/graph/parallel/iteration_macros.hpp>

#include <queue>
#include <algorithm> // for std::min, std::max
#include <iostream>
#include <atomic>
#include <tuple>
#include <climits>

template <typename OwnerMap, typename VertexComponentData>
struct vertex_component_owner {
  explicit vertex_component_owner(const OwnerMap& owner) : owner(owner) {}

  const OwnerMap& owner;
};

template <typename OwnerMap, typename VertexComponentData>
typename boost::property_traits<OwnerMap>::value_type
get(const vertex_component_owner<OwnerMap, VertexComponentData>& o, const VertexComponentData& data)
{ return get(o.owner, data.first); }

namespace boost {
namespace graph {
namespace distributed {

// The default priority queue
template<typename V, typename D, typename Compare>
struct default_cc_priority_queue {

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
struct cc_default_priority_queue_gen {
  template<typename V, typename D, typename Compare>
  struct queue {
    typedef default_cc_priority_queue<V, D, Compare> type;
  };
};


template<typename Graph,
	 typename ComponentMap,
	 typename IdDistribution,
         typename WorkStats,
	 typename PriorityQueueGenerator = cc_default_priority_queue_gen,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class data_driven_cc {
  // Type definitions for Vertex and Component
  // Vertex information is in graph_traits::vertex_descriptor
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

  typedef typename property_traits<ComponentMap>::value_type Component;
  // This is going to be our message
  typedef std::pair<Vertex, Component> vertex_component_data;
  typedef append_buffer<vertex_component_data, 10u> InitialBuffer;

  struct default_comparer {
    bool operator()(const vertex_component_data& vc1, const vertex_component_data& vc2) {
      if (vc1.second > vc2.second)
	return true;
      else if (vc1.second == vc2.second)
	return vc1.first > vc2.first;
      else
	return false;
    }
  };


  // Owner map type - To calculate which vertex belong to which node
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  // Create the default queue generator type
  typedef typename PriorityQueueGenerator::template queue<Vertex, Component, 
							  default_comparer>::type PriorityQueueType;

  // The handler, forward decleration
  struct vertex_component_handler;

  typedef typename MessageGenerator::template call_result<
    vertex_component_data,
    vertex_component_handler,
    vertex_component_owner<OwnerMap, vertex_component_data>,
    amplusplus::idempotent_combination_t<boost::parallel::minimum<Component>, Component> >::type
  RelaxMessage;

public:
  data_driven_cc(Graph& g,
                 amplusplus::transport &t,
		 int offs,
                 ComponentMap& component_map,
		 IdDistribution& dist,
                 WorkStats& stat,
                 MessageGenerator message_gen =
                 MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)),
		 unsigned int flushFreq = 2000)
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_component_data>(), 0)),
    graph(g),
    transport(t),
    core_offset(offs),
    component_map(component_map),
    id_distribution(dist),
    work_stats(stat),
    owner_map(get(vertex_owner, g)),
    relax_msg(message_gen,
	      transport,
	      vertex_component_owner<OwnerMap, vertex_component_data>(owner_map),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Component>(), std::numeric_limits<Component>::max())),
    flushFrequency(flushFreq)
  {
    initialize();
  }

  void operator() (int tid) { run(tid); }
  void run(int tid = 0);
  void populate_initial_buffer(int tid);

  template<typename SizeType>
  SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }

  time_type get_start_time() { return start_time; }
  time_type get_end_time() { return end_time; }

private:
  void initialize();
  void relax(const vertex_component_data&, const int);
  void handle_queue(const int, amplusplus::transport::end_epoch_request&);

  const int dummy_first_member_for_init_order; // Unused
  const Graph& graph;
  amplusplus::transport& transport;
  int core_offset;
  ComponentMap& component_map;
  IdDistribution& id_distribution;
  shared_ptr<InitialBuffer> buffer;
  WorkStats& work_stats;
  OwnerMap owner_map;
  unsigned int flushFrequency;
  RelaxMessage relax_msg;
  std::vector<PriorityQueueType> thread_buckets;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  unsigned int in_edge_factor;
  time_type start_time;
  time_type end_time;

}; // end of class

#define DD_CC_PARMS \
  typename Graph, typename ComponentMap, typename IdDistribution, typename WorkStats, typename PriorityQueueGenerator, typename MessageGenerator

#define DD_CC_TYPE \	
  data_driven_cc<Graph, ComponentMap, IdDistribution, WorkStats, PriorityQueueGenerator, MessageGenerator>

template<DD_CC_PARMS> void
DD_CC_TYPE::initialize() {
  relax_msg.set_handler(vertex_component_handler(*this));

  // Initialize components labels
  BGL_FORALL_VERTICES_T(v, graph, Graph) { 
    put(component_map, v, logical_id(v)); 
  }

  // Initial buffer
  shared_ptr<InitialBuffer> p(new InitialBuffer);
  buffer.swap(p);
}

template<DD_CC_PARMS>
struct DD_CC_TYPE::vertex_component_handler {
  vertex_component_handler() : self(NULL) {}
  vertex_component_handler(data_driven_cc& self) : self(&self) {}

  void operator() (const vertex_component_data& data) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->relax(data, tid);
  }

protected:
  data_driven_cc* self;
};


template<DD_CC_PARMS>
void
DD_CC_TYPE::populate_initial_buffer(int tid) {

  int nthreads = transport.get_nthreads();

  BGL_PARFORALL_VERTICES_T(u, graph, Graph, tid, nthreads) {
    Vertex min_neighbor = logical_id(u);
    std::set<Vertex> adjacencies1;
    BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
      Vertex v = target(e, graph);
      if (u != v) { // ignore self loops
	if (adjacencies1.insert(v).second) {
	  // check whether u is the minimum out of its neighbors
	  if (logical_id(v) < min_neighbor) {
	    min_neighbor = logical_id(v);
	    break;
	  }
	}
      }
    }

    if (min_neighbor == logical_id(u)) {
      buffer->push_back(vertex_component_data(u, logical_id(u)));
    }
  }
}


template<DD_CC_PARMS> void
DD_CC_TYPE::run(int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    int nthreads = transport.get_nthreads();

    if (0 == tid) {

      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));

      // Resize thread queue
      thread_buckets.resize(nthreads);

      for (int k = 0; k < nthreads; ++k) {
	// For each thread create a priority queue
	thread_buckets[k] = PriorityQueueType();
      }
    }

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

    populate_initial_buffer(tid);

    typename InitialBuffer::size_type buff_size;
    buff_size = buffer->size();


    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
    t_bar->wait();

    for (typename InitialBuffer::size_type i = tid ; 
	 i < buff_size ; i+= nthreads) {
      PriorityQueueType& queue = thread_buckets[tid];
      vertex_component_data& vd = (*buffer)[i];
      queue.put(vd);
    }

    // clear the initial buffer
    buffer->clear();

    { amplusplus::scoped_epoch epoch(transport); }

    // if we have work process queues. queue is populated in the relax_msg
    start_time = get_time();

    // Start the algorithm
    transport.begin_epoch();

    // Check whether we still have work to do
    amplusplus::transport::end_epoch_request request = transport.i_end_epoch();

    // increase activity count if queue is not empty
    PriorityQueueType& pqueue = thread_buckets[tid];
    if (!pqueue.is_empty())
      transport.increase_activity_count(1);

    handle_queue(tid, request);
    end_time = get_time();
  }

#ifdef PRINT_DEBUG
  std::cerr<<"Done with run\n";
#endif
}

template<DD_CC_PARMS>
void DD_CC_TYPE::handle_queue(const int tid, amplusplus::transport::end_epoch_request& request) {
  // get the queue
  PriorityQueueType& queue = thread_buckets[tid];
  int doFlushCounter = 0;

  while (true) {
    if (!queue.is_empty()) {
      // get the top vertex
      std::pair<Vertex, Component> vc = queue.pop();
      Vertex v = vc.first;
      Component v_component = vc.second;
      //      std::cout << "vertex : " << v << " component : " << v_component << std::endl;

      // remember if the queue was empty after pop
      bool was_empty = queue.is_empty();

      // get existing component
      // if component in component map is better that component we got from queue
      // then we dont need to do anything (i.e. another thread has updated component map with better
      // component), otherwise send new components for neighbours
      Component map_comp = get(component_map, v);

      bool haslowernbr = false;
      if (v_component == map_comp) {
	std::set<Vertex> adjacencies;
	BGL_FORALL_ADJ_T(v, u, graph, Graph) {
#ifdef PRINT_DEBUG
	  std::cout << "The v : " << v << " the u : " << u << std::endl;
#endif
	  if (logical_id(u) > v_component) { // ignore self-loops
	    adjacencies.insert(u);
	    //	    if (adjacencies.insert(u).second) {
	    // doFlushCounter++;
	    // relax_msg.send(vertex_component_data(u, v_component));
	    //}
	  } else if(logical_id(u) < v_component) { // v has a neighbor that is lower than v_component. Therefore, stop the search from v_component
	    haslowernbr = true;
	    break;
	  }
	}

	if (!haslowernbr) {
	  typename std::set<Vertex>::iterator ite = adjacencies.begin();
	  for(; ite != adjacencies.end(); ++ite) {
	    doFlushCounter++;
	    relax_msg.send(vertex_component_data((*ite), v_component));
	  }
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


template<DD_CC_PARMS> void
DD_CC_TYPE::relax(const vertex_component_data& vc, const int tid) {
  using boost::parallel::val_compare_and_swap;

  Vertex v = vc.first;
  Component new_component = vc.second;

  // get existing component
  // need last_old_dist to atomically update component map
  Component old_comp = get(component_map, v), last_old_comp;

#ifdef PRINT_DEBUG
  std::cout << "calling relax for vertex : " << v
	    << ", local id :" << logical_id(v)
	    << " new component : " << new_component 
	    << " old component : " << old_comp
	    << std::endl;

#endif

#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_edges(tid);
#endif

  // check whether we got a better component if we got better
  // component update the component map
  // credit : took from delta_stepping_shortest_paths.hpp
  while (new_component < old_comp) {
    last_old_comp = old_comp;
    old_comp = val_compare_and_swap(&component_map[v], old_comp, new_component);
    if (last_old_comp == old_comp) {

#ifdef PBGL2_PRINT_WORK_STATS
      if (last_old_comp < logical_id(v))
        work_stats.increment_invalidated(tid);
      else
        work_stats.increment_useful(tid);
#endif


      // we got a better component - insert it to work queue
      PriorityQueueType& queue = thread_buckets[tid];
      if (queue.is_empty()) {
	// If the queue was empty, increase activity count to notify AM++ that work is being held off.
#ifdef PRINT_DEBUG
	std::cout << "########### Increasing activity count by 1 ########## " << std::endl;
#endif
	transport.increase_activity_count(1);
      }
      queue.put(vc);
      return;
    }
  }

#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(tid);
#endif

}

} // end of distributed
} // end of graph
} // end of boost

#endif
