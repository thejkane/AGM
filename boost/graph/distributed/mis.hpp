// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Maximal Independent Set Algortihms================//
// This implements two algorithms: 1. FIX 
// MIS algorithm 2. FIX-PQ algorithm with thread level 
// ordering (use MIS_PRIORITY preprocessor macro)
//===========================================================//


#ifndef BOOST_GRAPH_MIS
#define BOOST_GRAPH_MIS

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/parallel/iteration_macros.hpp> // for all_reduce
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <algorithm> // for std::min, std::max
#include <boost/format.hpp>
#include <iostream>
#include <atomic>
#include "boost/tuple/tuple.hpp"
#include "thread_pq_def.hpp"
#include <boost/graph/distributed/owner_defs.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

//for profiling
#ifdef CRAYPAT
#include <pat_api.h>
#endif

#define MIS_UNFIX 0
#define MIS_FIX1 1
#define MIS_FIX0 2

typedef int state_t;

namespace boost { namespace graph { namespace distributed {

template<typename Graph, typename MISMap,
	 typename IdDistribution,
	 typename PriorityQueueGenerator = thread_priority_queue_gen,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class maximal_independent_set {

  typedef maximal_independent_set<Graph, MISMap, IdDistribution,
				  PriorityQueueGenerator, MessageGenerator> 
    self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;

  typedef append_buffer<Vertex, 10u> InitialBuffer;

  typedef unsigned long distance_t;

  // < Target, < < Source, Source State >, distance from min neighbor> >
  typedef std::pair<Vertex, std::pair<Vertex, std::pair<state_t, distance_t> > > work_item_t;

  static Vertex targetv(work_item_t wi) {
    return wi.first;
  }

  static Vertex sourcev(work_item_t wi) {
    return wi.second.first;
  }

  static state_t statev(work_item_t wi) {
    return wi.second.second.first;
  }

  static distance_t distancev(work_item_t wi) {
    return wi.second.second.second;
  }


  typedef typename std::map<Vertex, int > NeighborCountMap;

  struct default_comparer {
    bool operator()(const work_item_t& wi1, const work_item_t& wi2) {
      state_t s1 = statev(wi1);
      state_t s2 = statev(wi2);
      
      if (s1 > s2)
	return true;
#ifdef MIS_DISTANCE_ORDERING
      else if (s1 == s2) { // MIS_FIX1 vertices should be processed first
	distance_t d1 = distancev(wi1);
	distance_t d2 = distancev(wi2);

	if (d1 > d2) { 
	  return true;
	} else if (d1 == d2) {
	  Vertex v1 = targetv(wi1);
	  Vertex v2 = targetv(wi2);
	
	  if (v1 > v2)
	    return true;
	  else
	    return false;
	
	} else
	  return false;
      } else
	return false;
#else // distance ordering not defined
      else
	return false;
#endif

    }

    /*#ifdef 0
    bool operatorinverse()(const work_item_t& wi1, const work_item_t& wi2) const {
      state_t s1 = statev(wi1);
      state_t s2 = statev(wi2);
      
      if (s1 < s2)
	return true;
      else if (s1 == s2) {
	if (s1 == 0) {
	  distance_t d1 = distancev(wi1);
	  distance_t d2 = distancev(wi2);

	  if (d1 < d2) { // MIS_FIX1 vertices should be processed first
	    return true;
	  } else if (d1 == d2) {
	    Vertex v1 = targetv(wi1);
	    Vertex v2 = targetv(wi2);
	
	    if (v1 < v2)
	      return true;
	    else
	      return false;
	
	  } else
	    return false;
	} else {
	  Vertex v1 = targetv(wi1);
	  Vertex v2 = targetv(wi2);
	
	  if (v1 < v2)
	    return true;
	  else
	    return false;
	}
      } else
	return false;

    }
    #endif*/
  };


  typedef typename PriorityQueueGenerator::template queue<work_item_t, 
							  default_comparer>::type PriorityQueue_t;
  struct processing_function;

  struct minimum_pair_first
  {
    template<typename T>
    const T& operator()(const T& x, const T& y) const { return x.first < y.first ? x : y; }

    template<typename F>
    struct result {
      typedef typename boost::function_traits<F>::arg1_type type;
    };
  };


  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;
  
public:
  maximal_independent_set(Graph& g,
			  MISMap& mismap, 
			  amplusplus::transport &t,
			  const IdDistribution& idd,
			  int freq,
			  int offset,
			  MessageGenerator message_gen = 
			  MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
      g(g), transport(t), 
      id_distribution(idd),
      mis(mismap),
      owner(get(vertex_owner, g)), 
      core_offset(offset),
    relax_msg(message_gen, transport, owner_from_pair<OwnerMap, work_item_t>(owner),
	      amplusplus::idempotent_combination(minimum_pair_first())),
      flushFrequency(freq), pq(t.get_nthreads()), counter(0)
  {
    initialize();
  }

  void operator() (int tid) { 

    run(tid); 

  }

  void run(int tid = 0);
  time_type get_start_time() {
    return start_time;
  }


  void add_to_pq(const work_item_t& data, int tid) {
    pq.push(data, tid);
  }

#ifdef MIS_STATS
  void print_stats() {
    unsigned long skipped = 0;
    unsigned long called = 0;

    std::vector<unsigned long>::iterator ite = skipped_handlers.begin();
    for (; ite != skipped_handlers.end(); ++ite) {
      skipped += (*ite);
    }

    ite = called_handlers.begin();
    for (; ite != called_handlers.end(); ++ite) {
      called += (*ite);
    }
    
    std::cout << "Skiped : " << skipped << ", Called : " << called
	      << std::endl;
  }
#endif

#ifdef MIS_PRIORITY
  void print_pq_sizes() {
#ifdef MIS_STATS
    int i = 0;
    std::vector<size_t>::iterator ite = pq_sizes.begin();
    for(; ite != pq_sizes.end(); ++ite) {
      std::cout << "thread : " << i << ", pq size : " << (*ite)
		<< std::endl;
      ++i;
    }
#endif
  }
#endif


private:
  void initialize();
  void process(const work_item_t& data, int tid);
  void handle_queue(const int tid, amplusplus::transport::end_epoch_request& request);
  void handle_fix0_wi(const work_item_t& data, int tid, int& flushcnt);

  template<typename SizeType>
  inline SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }

  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   state_t state,
			   distance_t distance) {
    // typedef std::pair<Vertex, std::pair<Vertex, std::pair<state_t, distance_t> > > work_item_t;
    work_item_t wi(d, std::make_pair(s, std::make_pair(state, distance)));
    return wi;
  }

  const int dummy_first_member_for_init_order;
  const Graph& g;
  amplusplus::transport& transport;
  const IdDistribution& id_distribution;
  MISMap& mis;
  OwnerMap owner;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  RelaxMessage relax_msg;
  time_type start_time;
  int flushFrequency;
  PriorityQueue_t pq;
  std::atomic_ullong counter;
  NeighborCountMap lower_neighbors;
  NeighborCountMap lower_fixed_neighbors;
  shared_ptr<InitialBuffer> buffer;
  int core_offset;
#ifdef MIS_STATS
  std::vector<unsigned long> skipped_handlers;
  std::vector<unsigned long> called_handlers;
#endif
#ifdef MIS_PRIORITY
#ifdef MIS_STATS
  std::vector<size_t> pq_sizes;
#endif
#endif

};

#define MIS_PARMS                                   \
      typename Graph, typename MISMap, typename IdDistribution, typename PriorityQueueGenerator, typename MessageGenerator

#define MIS_TYPE                                    \
      maximal_independent_set<Graph, MISMap, IdDistribution, PriorityQueueGenerator, MessageGenerator>


template<MIS_PARMS>
void
MIS_TYPE::initialize() {

#ifdef MIS_PRIORITY
  std::cout << "Running with MIS priority ....." << std::endl;
#ifdef MIS_STATS
  pq_sizes.resize(transport.get_nthreads(), 0);
#endif
#else
  std::cout << "Running with *out* MIS priority ....." << std::endl;
#endif

  relax_msg.set_handler(processing_function(*this));

#ifdef MIS_STATS
  skipped_handlers.resize(transport.get_nthreads(), 0);
  called_handlers.resize(transport.get_nthreads(), 0);
#endif

  //
  shared_ptr<InitialBuffer> p(new InitialBuffer);
  buffer.swap(p);

  BGL_FORALL_VERTICES_T(u, g, Graph) {
    mis[u] = MIS_UNFIX;
    // lower_neighbors are not initialized
    std::set<Vertex> locadjacencies;
    int count = 0;
    BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
      Vertex v = target(e, g);
      if (logical_id(v) < logical_id(u)) {
	if (locadjacencies.insert(v).second) { //to avoid parallel edges
	  ++count;
	}
      }
    }

    if (count == 0) {
      buffer->push_back(u);
    }

    lower_neighbors.insert(std::make_pair(u, count));
    lower_fixed_neighbors.insert(std::make_pair(u, 0));
  }

  

}


template<MIS_PARMS>
void
MIS_TYPE::handle_fix0_wi(const work_item_t& data, int tid, int& flushcnt) {
  Vertex dest_vertex = targetv(data);
  Vertex source_vertex = sourcev(data);
  state_t source_state = statev(data);
  distance_t dist = distancev(data);

#ifdef MIS_PRIORITY
  if (mis[dest_vertex] != MIS_UNFIX) {
#ifdef PRINT_DEBUG
    assert(mis[dest_vertex] == MIS_FIX1 ||
	   mis[dest_vertex] == MIS_FIX0);

    if (mis[dest_vertex] == MIS_FIX1) {
      if (logical_id(source_vertex) < logical_id(dest_vertex)) {
	std::cout << "To dest_vertex : " << dest_vertex << " to be fixed we must have checked all its lower neighbors (" 
		  << source_vertex << ")" <<
	  "lfn=" << lower_fixed_neighbors[dest_vertex] << " ln=" << lower_neighbors[dest_vertex] << 
	  std::endl;
	assert(false);
      }
    }
#endif

#ifdef MIS_STATS
    ++skipped_handlers[tid];
#endif
    return;
  }
#endif

  // case 3:
  // Source vertex is not in mis and fixed =>
  // If source vertex is less than destination vertex, then add
  // 1 to lower_fixed_neighbors
  // Also check, all neighbors less than current vertex are in
  // lower_fixed_neighbors, if so promote dest_vertex to be in MIS and fix (MIS_FIX1).
  // Inform neigbors about state changes.
  assert(source_state == MIS_FIX0);

#ifdef MIS_STATS
  ++called_handlers[tid];
#endif

  if (logical_id(source_vertex) < logical_id(dest_vertex)) {

#ifdef PRINT_DEBUG
    assert(lower_neighbors[dest_vertex] > 0);
    if (lower_fixed_neighbors[dest_vertex] >= lower_neighbors[dest_vertex]) {
      std::cout << "lfn=" << lower_fixed_neighbors[dest_vertex] << " ln=" << lower_neighbors[dest_vertex]
		<< std::endl;
    }
#endif

    // if all lower neigbors are fixed dest vertex must also be fixed
    assert(lower_fixed_neighbors[dest_vertex] < lower_neighbors[dest_vertex]);

    if (__atomic_add_fetch(&lower_fixed_neighbors[dest_vertex], 1, __ATOMIC_SEQ_CST) == lower_neighbors[dest_vertex]) {
      // all the lower neigbors are fixed
      assert(lower_fixed_neighbors[dest_vertex] == lower_neighbors[dest_vertex]);
      mis[dest_vertex] = MIS_FIX1;
      // vertex is fixed. We do not need to maintain lower_fixed_neighbors
      // for this vertex any more.
      // source not in mis -- just inform every neighbor the current state
      std::set<Vertex> adjacencies;
      BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
	Vertex u = target(e, g);
	if (logical_id(u) > logical_id(dest_vertex)) { // we know that all lower neighbors are fixed, therefore we only need to send messages to higher neighbors
	  if (adjacencies.insert(u).second) {
#ifdef MIS_PRIORITY
	    //	    if (fix0vertices[tid].insert(u).second) {
	      ++flushcnt;
#endif
	      work_item_t wi = construct_wi(u, dest_vertex, mis[dest_vertex],
					    (dist + 1));
#ifndef MIS_PRIORITY
	      //	      transport.increase_activity_count(1);
#endif
	      relax_msg.send(wi);
#ifdef MIS_PRIORITY
	      // }
#endif
	  }
	}
      }
    }
  }
}


template<MIS_PARMS>
void
MIS_TYPE::process(const work_item_t& data, int tid) {
  //   typedef std::pair<Vertex, std::pair<Vertex, state_t > > work_item_t;

  /*#ifdef CRAYPAT
  if (PAT_region_begin ( 2, "misprocess" ) == PAT_API_FAIL) {
    std::cout << "PAT begin failed ! " << std::endl;
    assert(false);
  }
  #endif*/
  

  Vertex dest_vertex = targetv(data);
  Vertex source_vertex = sourcev(data);
  state_t source_state = statev(data);
  distance_t dist = distancev(data);


#ifdef PRINT_DEBUG
  std::cout << "dv : " << dest_vertex << " sv : " << source_vertex << " s_state : " << source_state
	    << " d_state :" << mis[dest_vertex] << " distance : " << dist 
	    << std::endl;
#endif

  // case 1:
  // is destination vertex fixed ? then just ignore
  if (mis[dest_vertex] != MIS_UNFIX) {
    assert(mis[dest_vertex] == MIS_FIX1 ||
	   mis[dest_vertex] == MIS_FIX0);
#ifndef MIS_PRIORITY
    //    transport.decrease_activity_count(1);
#endif

    /*#ifdef CRAYPAT
    if (PAT_region_end(2) == PAT_API_FAIL) {
     std::cout << "PAT end failed ! " << std::endl;
     assert(false);
    }
    #endif*/

#ifdef MIS_STATS
    ++skipped_handlers[tid];
#endif
    return;
  }

  // case 2:
  // Source vertex is in mis and fixed =>
  // Destination vertex state must be transferred to MIS_FIX0.
  if (source_state == MIS_FIX1) {
#ifdef MIS_STATS
    ++called_handlers[tid];
#endif
    //mis[dest_vertex] = MIS_FIX0; 
    // source in mis, therefore dest must be out of mis
    int expected = MIS_UNFIX;
    if (__atomic_compare_exchange_n(&mis[dest_vertex], &expected, MIS_FIX0, false /*no weak*/, 
				    __ATOMIC_SEQ_CST, __ATOMIC_RELAXED)) {

      assert(mis[dest_vertex] == MIS_FIX0);
      // inform all neighbors
      std::set<Vertex> adjacencies;
      BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
	Vertex u = target(e, g);
	if (logical_id(u) > logical_id(dest_vertex)) {
	  if (adjacencies.insert(u).second) {
	    work_item_t wi = construct_wi(u, dest_vertex, mis[dest_vertex],
					  (dist + 1));
	      
#ifdef MIS_PRIORITY
	    relax_msg.send(wi);
#else
	    //	    transport.increase_activity_count(1);
	    relax_msg.send(wi);
#endif
	  }
	}
      }
      
    }
  } else {
#ifdef MIS_PRIORITY

    if (pq.empty(tid)) {
      transport.increase_activity_count(1);
    }
    pq.put(data, tid);

#else
    int ignore = 0;
    handle_fix0_wi(data, tid, ignore);
#endif
  } // end of else

#ifndef MIS_PRIORITY
  // done processing the work item
  //  transport.decrease_activity_count(1);
#endif

  /*#ifdef CRAYPAT
  if (PAT_region_end(2) == PAT_API_FAIL) {
    std::cout << "PAT end failed ! " << std::endl;
    assert(false);
  }
  #endif*/


}

template<MIS_PARMS>
void
MIS_TYPE::run(int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {

#ifdef CRAYPAT
      if (PAT_region_begin ( 1, "misrun" ) == PAT_API_FAIL) {
	std::cout << "PAT begin failed ! " << std::endl;
	assert(false);
      }
#endif

    int nthreads = transport.get_nthreads();
    if (0 == tid) {
      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    }

    { amplusplus::scoped_epoch epoch(transport); }

    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
    t_bar->wait();

    // if two processes are running on the same node, core_offset
    // is important to achieve thread affinity
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

    typename InitialBuffer::size_type buff_size;
    buff_size = buffer->size();
    for (typename InitialBuffer::size_type i = tid ; 
	 i < buff_size ; i+= nthreads) {
      Vertex u = (*buffer)[i];
      mis[u] = MIS_FIX1;
      std::set<Vertex> adjacencies;
      BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
	Vertex v = target(e, g);
	if (u != v) { // ignore self loops
	  if (adjacencies.insert(v).second) {
	    work_item_t wi = construct_wi(v, u, mis[u], 0);
	    relax_msg.send(wi);
	  }
	}
      }

    }
#ifdef PRINT_DEBUG
    // for information, can remove
    if (tid == 0)
      std::cout << "Number of minimum_neighbor_vertices : " << counter << std::endl;
#endif

    amplusplus::transport::end_epoch_request request = transport.i_end_epoch();
    // if we have work process queues. queue is populated in the relax_msg
    handle_queue(tid, request);

#ifdef CRAYPAT
    if (PAT_region_end(1) == PAT_API_FAIL) {
      std::cout << "PAT end failed ! " << std::endl;
      assert(false);
    }
#endif



  }
}

template<MIS_PARMS>
void 
MIS_TYPE::handle_queue(const int tid, amplusplus::transport::end_epoch_request& request) {
  int doFlushCounter = 0;
  while(true) {
#ifdef MIS_PRIORITY
#ifdef MIS_STATS
    if (pq.size(tid) > pq_sizes[tid]) {
      pq_sizes[tid] = pq.size(tid);
    }
#endif

    work_item_t wi;
    while(pq.pop(wi, tid)) {
      // to make sure pq is ordering properly TODO remove once verified
#ifdef PRINT_DEBUG
      std::cout << "target vertex : " << targetv(wi) << " source : " << sourcev(wi) << std::endl;
#endif

      bool e = pq.empty(tid);
      handle_fix0_wi(wi, tid, doFlushCounter);

      if (e)
	transport.decrease_activity_count(1);

      if(doFlushCounter == flushFrequency) {
	doFlushCounter = 0;
	if (request.test()) {
	  assert(pq.empty(tid));
	  break;
	}
      }
    }

    if (request.test())
      break;
#else
    if (request.test())
      break;
#endif
  }
}

template<MIS_PARMS>
struct MIS_TYPE::
processing_function {
  
  processing_function() : self(NULL) {}
  processing_function(maximal_independent_set& self) : self(&self) {}
  
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
#ifdef PRINT_DEBUG
    std::cout << "Handler called in tid : " << tid << std::endl;
#endif
    self->process(data, tid);
  }

protected:
  maximal_independent_set* self;
};

}}}
#endif
