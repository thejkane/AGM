// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef DATA_DRIVEN_DIFF_CONNECTED_COMPONENTS
#define DATA_DRIVEN_DIFF_CONNECTED_COMPONENTS
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
	 typename ComponentMap,
	 typename PriorityQueueGenerator = cc_default_priority_queue_gen,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class data_diff_driven_cc {
  // Type definitions for Vertex and Component
  // Vertex information is in graph_traits::vertex_descriptor
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

  typedef typename property_traits<ComponentMap>::value_type Component;
  // This is going to be our message
  typedef std::pair<Vertex, Component> vertex_component_data;

  struct default_comparer {
    bool operator()(const vertex_component_data& vc1, const vertex_component_data& vc2) {
      if (vc1.second > vc2.second)
	return true;
      else if (vc1.second == vc2.second)
	return vc1.first > vc2.first;
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
  data_diff_driven_cc(Graph& g,
		      ComponentMap component_map,
      		      amplusplus::transport &t,
		      MessageGenerator message_gen =
		        MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)),
		      MessageGenerator priority_message_gen =
		        MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12), 1),
		      unsigned int flushFreq = 2000,
		      unsigned int eagerLimit = 1000)
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_component_data>(), 0)),
    graph(g),
    transport(t),
    component_map(component_map),
    owner_map(get(vertex_owner, g)),
    relax_msg(message_gen,
	      transport,
	      vertex_component_owner<OwnerMap, vertex_component_data>(owner_map),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Component>(), std::numeric_limits<Component>::max())),
    priority_relax_msg(priority_message_gen,
	      transport,
	      vertex_component_owner<OwnerMap, vertex_component_data>(owner_map),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Component>(), std::numeric_limits<Component>::max())),
    flushFrequency(flushFreq),
    eagerLimit(eagerLimit),
    counter(0)
  {
    initialize();
  }

  void operator() (int tid) { run(tid); }
  void run(int tid = 0);
  unsigned int get_counter() {
    return counter;
  }


private:
  void initialize();
  void relax(const vertex_component_data&, const int);
  void handle_queue(const int, amplusplus::transport::end_epoch_request&);

  const int dummy_first_member_for_init_order; // Unused
  const Graph& graph;
  amplusplus::transport& transport;
  ComponentMap component_map;
  OwnerMap owner_map;
  unsigned int flushFrequency;
  unsigned int eagerLimit;
  RelaxMessage relax_msg;
  RelaxMessage priority_relax_msg;
  std::vector<PriorityQueueType> thread_buckets;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  unsigned int counter;

}; // end of class

#define DD_CC_PARMS \
  typename Graph, typename ComponentMap, typename PriorityQueueGenerator, typename MessageGenerator

#define DD_CC_TYPE \	
  data_diff_driven_cc<Graph, ComponentMap, PriorityQueueGenerator, MessageGenerator>

template<DD_CC_PARMS> void
DD_CC_TYPE::initialize() {
  relax_msg.set_handler(vertex_component_handler(*this));
  priority_relax_msg.set_handler(vertex_component_handler(*this));
}

template<DD_CC_PARMS>
struct DD_CC_TYPE::vertex_component_handler {
  vertex_component_handler() : self(NULL) {}
  vertex_component_handler(data_diff_driven_cc& self) : self(&self) {}

  void operator() (const vertex_component_data& data) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->relax(data, tid);
  }

protected:
  data_diff_driven_cc* self;
};

template<DD_CC_PARMS> void
DD_CC_TYPE::run(int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    if (0 == tid) {
      int nthreads = transport.get_nthreads();

      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));

      // Resize thread queue
      thread_buckets.resize(nthreads);

      for (int k = 0; k < nthreads; ++k) {
	// For each thread create a priority queue
	thread_buckets[k] = PriorityQueueType();
      }
    }

    // TODO delta had following. Not sure we also need this - ask
    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized
    { amplusplus::scoped_epoch epoch(transport); }

    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
    t_bar->wait();

    // Start the algorithm
    transport.begin_epoch();

    /*    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    
    std::pair<vertex_iter, vertex_iter> vp;
    for (vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
      if (get(owner_map, *vp.first) == transport.rank()) {
	relax(vertex_component_data(*vp.first, *vp.first), tid);
      }
      }*/

    BGL_FORALL_VERTICES_T(v, graph, Graph) {
      if (get(owner_map, v) == transport.rank()) {
	relax(vertex_component_data(v, v), tid);
      }
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

template<DD_CC_PARMS>
void DD_CC_TYPE::handle_queue(const int tid, amplusplus::transport::end_epoch_request& request) {
  // get the queue
  PriorityQueueType& queue = thread_buckets[tid];

  while (true) {
    if (!queue.is_empty()) {
      // get the top vertex
      std::pair<Vertex, Component> vc = queue.pop();
      Vertex v = vc.first;
      Component v_component = vc.second;

      // remember if the queue was empty after pop
      bool was_empty = queue.is_empty();

      // get existing component
      // if component in component map is better that component we got from queue
      // then we dont need to do anything (i.e. another thread has updated component map with better
      // component), otherwise send new components for neighbours
      Component map_comp = get(component_map, v);

      if (v_component == map_comp) {
	BGL_FORALL_ADJ_T(v, u, graph, Graph) {
#ifdef PRINT_DEBUG
	  std::cout << "The v : " << v << " the u : " << u << std::endl;
#endif
	  relax_msg.send(vertex_component_data(u, v_component));
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

      if(request.test()) {
#ifdef PRINT_DEBUG
	assert(queue.is_empty());
#endif
	break;
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
	    << " new component : " << new_component 
	    << " old component : " << old_comp
	    << std::endl;
#endif

  if (old_comp == std::numeric_limits<Vertex>::max() || old_comp == v) {
      while (new_component < old_comp) {
	last_old_comp = old_comp;
	old_comp = val_compare_and_swap(&component_map[v], old_comp, new_component);
	if (last_old_comp == old_comp) {
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
    } else if(new_component < old_comp) {
      counter ++;
      // RECORD RE-WRITE (old_comp, new_comp)
    }
}

} // end of distributed
} // end of graph
} // end of boost

#endif
