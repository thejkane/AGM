// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Connected Components Algortihm================//
// Chaotic version of connected components algorithm.
//===========================================================//

#ifndef DATA_DRIVEN_CONNECTED_COMPONENTS_CHAOTIC
#define DATA_DRIVEN_CONNECTED_COMPONENTS_CHAOTIC
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
	 typename IdDistribution,
         typename WorkStats,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class chaotic_cc {
  // Type definitions for Vertex and Component
  // Vertex information is in graph_traits::vertex_descriptor
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type VerticesSize;
  typedef typename property_traits<ComponentMap>::value_type Component;
  // This is going to be our message
  typedef std::pair<Vertex, Component> vertex_component_data;

  typedef append_buffer<vertex_component_data, 10u> InitialBuffer;

  // Owner map type - To calculate which vertex belong to which node
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  // The handler, forward decleration
  struct vertex_component_handler;

  typedef typename MessageGenerator::template call_result<
    vertex_component_data,
    vertex_component_handler,
    vertex_component_owner<OwnerMap, vertex_component_data>,
    amplusplus::idempotent_combination_t<boost::parallel::minimum<Component>, Component> >::type
  RelaxMessage;

public:
  chaotic_cc(Graph& g,
	     amplusplus::transport &t,
	     int offs,
	     ComponentMap& component_map,
	     IdDistribution& dist,
	     WorkStats& stat,
	     MessageGenerator message_gen =
	     MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
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
	      amplusplus::idempotent_combination(boost::parallel::minimum<Component>(), std::numeric_limits<Component>::max()))
  {
    initialize();
  }

  void operator() (int tid) { run(tid); }
  void run(int tid = 0);
  void populate_initial_buffer(int tid);

  time_type get_start_time() { return start_time; }
  time_type get_end_time() { return end_time; }


  template<typename SizeType>
  SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }


private:
  void initialize();
  void relax(const vertex_component_data&, const int);

  const int dummy_first_member_for_init_order; // Unused
  const Graph& graph;
  amplusplus::transport& transport;
  int core_offset;
  ComponentMap component_map;
  IdDistribution& id_distribution;
  shared_ptr<InitialBuffer> buffer;
  WorkStats& work_stats;
  OwnerMap owner_map;
  RelaxMessage relax_msg;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  time_type start_time;
  time_type end_time;

}; // end of class

#define CHAOTIC_CC_PARMS \
  typename Graph, typename ComponentMap, typename IdDistribution, typename WorkStats, typename MessageGenerator

#define CHAOTIC_CC_TYPE \	
  chaotic_cc<Graph, ComponentMap, IdDistribution, WorkStats, MessageGenerator>

template<CHAOTIC_CC_PARMS> void
CHAOTIC_CC_TYPE::initialize() {
  relax_msg.set_handler(vertex_component_handler(*this));
  // Initialize components labels
  BGL_FORALL_VERTICES_T(v, graph, Graph) { 
    put(component_map, v, logical_id(v)); 
  }

  // Initial buffer
  shared_ptr<InitialBuffer> p(new InitialBuffer);
  buffer.swap(p);
}

template<CHAOTIC_CC_PARMS>
struct CHAOTIC_CC_TYPE::vertex_component_handler {
  vertex_component_handler() : self(NULL) {}
  vertex_component_handler(chaotic_cc& self) : self(&self) {}

  void operator() (const vertex_component_data& data) const {
    const int tid = amplusplus::detail::get_thread_id();
    self->relax(data, tid);
  }

protected:
  chaotic_cc* self;
};

template<CHAOTIC_CC_PARMS>
void
CHAOTIC_CC_TYPE::populate_initial_buffer(int tid) {

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
      std::set<Vertex> adjacencies;
      BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
	Vertex v = target(e, graph);
	if (v != u) {
	  if (adjacencies.insert(v).second) {
	    buffer->push_back(vertex_component_data(v, logical_id(u)));
	  }
	}
      }
    }
  }
}


template<CHAOTIC_CC_PARMS> void
CHAOTIC_CC_TYPE::run(int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    int nthreads = transport.get_nthreads();

    if (0 == tid) {
      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    }

    // TODO delta had following. Not sure we also need this - ask
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

    validate_thread_core_relation();

    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
    t_bar->wait();

    populate_initial_buffer(tid);

    { amplusplus::scoped_epoch epoch(transport); }

    typename InitialBuffer::size_type buff_size;
    buff_size = buffer->size();

    // start timer
    start_time = get_time();

    { 
      amplusplus::scoped_epoch epoch(transport);
      for (typename InitialBuffer::size_type i = tid ; 
	   i < buff_size ; i+= nthreads) {
	vertex_component_data& vd = (*buffer)[i];
	relax_msg.send(vd);
      }
    }

    // clear the initial buffer
    buffer->clear();

    end_time = get_time();
  }

#ifdef PRINT_DEBUG
  std::cerr<<"Done with run\n";
#endif
}


template<CHAOTIC_CC_PARMS> void
CHAOTIC_CC_TYPE::relax(const vertex_component_data& vc, const int tid) {
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

      bool haslowernbr = false;
      std::set<Vertex> adjacencies;
      BGL_FORALL_ADJ_T(v, u, graph, Graph) {
	Vertex lu = logical_id(u);
#ifdef PRINT_DEBUG
	std::cout << "The v : " << v << " the u : " << u << std::endl;
#endif
	if (lu > new_component) {
	  if (v != u) { // ignore self-loops
	    adjacencies.insert(u);
	  }
	} else if (lu < new_component) {
	  // has a lower neighbor
	  haslowernbr = true;
	  break;
	}
      }

      if (!haslowernbr) {
	typename std::set<Vertex>::iterator ite = adjacencies.begin();
	for(; ite != adjacencies.end(); ++ite) {
	  relax_msg.send(vertex_component_data((*ite), new_component));
	}
      }
    }

    return;
  }

#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(tid);
#endif

}

} // end of distributed
} // end of graph
} // end of boost

#endif
