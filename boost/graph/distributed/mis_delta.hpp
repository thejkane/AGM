// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Maximal Independent Set Algortihms================//
// AGM variation of MIS algorithm.
//===========================================================//

#ifndef BOOST_GRAPH_MIS_DELTA
#define BOOST_GRAPH_MIS_DELTA

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
#include <boost/graph/distributed/owner_defs.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

#define MIS_UNFIX 0
#define MIS_FIX1 1
#define MIS_FIX0 2

typedef int state_t;

namespace boost { namespace graph { namespace distributed {

template<typename Graph, typename MISMap,
	 typename IdDistribution,
	 typename Bucket = append_buffer<typename graph_traits<Graph>::vertex_descriptor, 10u>,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class maximal_independent_set_delta {

  typedef maximal_independent_set_delta<Graph, MISMap, IdDistribution, Bucket,
				  MessageGenerator> self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;

  typedef typename std::vector<Bucket*>::size_type BucketIndex;

  // < Target, < Source, Source State > >
  typedef std::pair<Vertex, std::pair<Vertex, state_t> > work_item_t;

  static Vertex targetv(const work_item_t& wi) {
    return wi.first;
  }

  static Vertex sourcev(const work_item_t& wi) {
    return wi.second.first;
  }

  static state_t statev(const work_item_t& wi) {
    return wi.second.second;
  }

  typedef typename std::map<Vertex, int > NeighborCountMap;

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
  maximal_independent_set_delta(Graph& g,
				MISMap& mismap, 
				amplusplus::transport &t,
				const IdDistribution& idd,
				int freq,
				int offset,
				MessageGenerator message_gen = 
				MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
      g(g), mis(mismap), transport(t), 
      id_distribution(idd),
      owner(get(vertex_owner, g)), 
      core_offset(offset),
    relax_msg(message_gen, transport, owner_from_pair<OwnerMap, work_item_t>(owner),
	      amplusplus::idempotent_combination(minimum_pair_first())),
      flushFrequency(freq), counter(0)
  {
    initialize();
  }

  void operator() (int tid) { run(tid); }
  void run(int tid = 0);

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


  time_type get_start_time() {
    return start_time;
  }


private:
  void initialize();
  void process(const work_item_t& data, int tid);
  void handle_buckets(const int tid, int nthreads);
  void handle_fix0_wi(const work_item_t& data, int tid);

  template<typename SizeType>
  inline SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }


  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   state_t state) {
    // typedef std::pair<Vertex, std::pair<Vertex, state_t > > work_item_t;
    work_item_t wi(d, std::make_pair(s, state));
    return wi;
  }

  const int dummy_first_member_for_init_order;
  const Graph& g;
  MISMap& mis;
  amplusplus::transport& transport;
  const IdDistribution& id_distribution;
  OwnerMap owner;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  RelaxMessage relax_msg;
  time_type start_time;
  int flushFrequency;
  std::atomic_ullong counter;
  NeighborCountMap lower_neighbors;
  NeighborCountMap lower_fixed_neighbors;
  std::vector<shared_ptr<Bucket> > buckets;
  int core_offset;
#ifdef MIS_STATS
  std::vector<unsigned long> skipped_handlers;
  std::vector<unsigned long> called_handlers;
#endif

};

#define MIS_PARMS_DELTA                                   \
      typename Graph, typename MISMap, typename IdDistribution, typename Bucket, typename MessageGenerator

#define MIS_TYPE_DELTA                                    \
      maximal_independent_set_delta<Graph, MISMap, IdDistribution, Bucket, MessageGenerator>


#define FIX1_INDEX 0
#define FIX0_INDEX 1

template<MIS_PARMS_DELTA>
void
MIS_TYPE_DELTA::initialize() {

  relax_msg.set_handler(processing_function(*this));
  // only two buckets
  // bucket 0 - fix1 vertices
  // bucket 1 - fix0 vertices
  buckets.resize(2);
  for (BucketIndex i = 0 ; i < buckets.size() ; ++i) {
    shared_ptr<Bucket> p(new Bucket);
    buckets[i].swap(p);
  }

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

    lower_neighbors.insert(std::make_pair(u, count));
    lower_fixed_neighbors.insert(std::make_pair(u, 0));
  }

#ifdef MIS_STATS
  skipped_handlers.resize(transport.get_nthreads(), 0);
  called_handlers.resize(transport.get_nthreads(), 0);
#endif


}

template<MIS_PARMS_DELTA>
void
MIS_TYPE_DELTA::handle_fix0_wi(const work_item_t& data, int tid) {
  Vertex dest_vertex = targetv(data);
  Vertex source_vertex = sourcev(data);
  state_t source_state = statev(data);

  if (mis[dest_vertex] != MIS_UNFIX) {
    assert(mis[dest_vertex] == MIS_FIX1 ||
	   mis[dest_vertex] == MIS_FIX0);

#ifdef MIS_STATS
    ++skipped_handlers[tid];
#endif

    return;
  }


  // case 3:
  // Source vertex is not in mis and fixed =>
  // If source vertex is less than destination vertex, then add
  // 1 to lower_fixed_neighbors
  // Also check, all neighbors less than current vertex are in
  // lower_fixed_neighbors, if so promote dest_vertex to be in MIS and fix (MIS_FIX1).
  // Inform neigbors about state changes.
  assert(source_state == MIS_FIX0);

  if (logical_id(source_vertex) < logical_id(dest_vertex)) {

#ifdef MIS_STATS
  ++called_handlers[tid];
#endif


#ifdef PRINT_DEBUG
    assert(lower_neighbors[dest_vertex] > 0);
    if (lower_fixed_neighbors[dest_vertex] >= lower_neighbors[dest_vertex]) {
      std::cout << "lfn=" << lower_fixed_neighbors[dest_vertex] << " ln=" << lower_neighbors[dest_vertex]
		<< std::endl;
    }

    std::cout << "[BEFORE] lfn=" << lower_fixed_neighbors[dest_vertex] << " ln=" << lower_neighbors[dest_vertex]
		<< std::endl;
#endif

    // if all lower neigbors are fixed dest vertex must also be fixed
    assert(lower_fixed_neighbors[dest_vertex] < lower_neighbors[dest_vertex]);

    if (__atomic_add_fetch(&lower_fixed_neighbors[dest_vertex], 1, __ATOMIC_SEQ_CST) 
	== lower_neighbors[dest_vertex]) {
#ifdef PRINT_DEBUG
      // all the lower neigbors are fixed
      std::cout << "============ Inside If Condition =============" << std::endl;
#endif
      assert(lower_fixed_neighbors[dest_vertex] == lower_neighbors[dest_vertex]);
      mis[dest_vertex] = MIS_FIX1;
      buckets[FIX1_INDEX]->push_back(dest_vertex);
    }

#ifdef PRINT_DEBUG
    std::cout << "[AFTER] lfn=" << lower_fixed_neighbors[dest_vertex] << " ln=" << lower_neighbors[dest_vertex]
		<< std::endl;
#endif

  }
}


template<MIS_PARMS_DELTA>
void
MIS_TYPE_DELTA::process(const work_item_t& data, int tid) {
  
  Vertex dest_vertex = targetv(data);
  Vertex source_vertex = sourcev(data);
  state_t source_state = statev(data);


#ifdef PRINT_DEBUG
  std::cout << "dv : " << dest_vertex << " sv : " << source_vertex << " s_state : " << source_state
	    << " d_state :" << mis[dest_vertex]
	    << std::endl;
#endif

  // case 1:
  // is destination vertex fixed ? then just ignore -- Ideally we should not hit this!
  if (mis[dest_vertex] != MIS_UNFIX) {
    assert(mis[dest_vertex] == MIS_FIX1 ||
	   mis[dest_vertex] == MIS_FIX0);

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
      buckets[FIX0_INDEX]->push_back(dest_vertex);      
    }
  } else {
    handle_fix0_wi(data, tid);
  } // end of else
}

template<MIS_PARMS_DELTA>
void
MIS_TYPE_DELTA::run(int tid) {
 
  //  print_tid_coreid();

  AMPLUSPLUS_WITH_THREAD_ID(tid) {

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
    
    { 
      amplusplus::scoped_epoch epoch(transport); 

      BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
	if (get(owner, u) == transport.rank()) {

	  Vertex min_neighbor = u;
	  std::set<Vertex> adjacencies1;
	  BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
	    Vertex v = target(e, g);
	    if (u != v) { // ignore self loops
	      if (adjacencies1.insert(v).second) {
		// check whether u is the minimum out of its neighbors
		if (logical_id(v) < logical_id(min_neighbor)) {
		  min_neighbor = v;
		  break;
		}
	      }
	    }
	  }

	  if (logical_id(min_neighbor) == logical_id(u)) {
#ifdef PRINT_DEBUG
	    counter++; // for information. can remove
#endif
	    mis[u] = MIS_FIX1;
	    std::set<Vertex> adjacencies;
	    BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
	      Vertex v = target(e, g);
	      if (u != v) { // ignore self loops
		if (adjacencies.insert(v).second) {
		  work_item_t wi = construct_wi(v, u, mis[u]);
		  relax_msg.send(wi);
		}
	      }
	    }
	  }
	}
      }
    } // end of scoped epoch
    
#ifdef PRINT_DEBUG
    // for information, can remove
    if (tid == 0)
      std::cout << "Number of minimum_neighbor_vertices : " << counter << std::endl;
#endif

    // if we have work process queues. queue is populated in the relax_msg
    handle_buckets(tid, nthreads);
  }

  //{ amplusplus::scoped_epoch epoch(transport); }
  //std::cout << "done .." << std::endl;
  //print_tid_coreid();
}

template<MIS_PARMS_DELTA>
void 
MIS_TYPE_DELTA::handle_buckets(const int tid, int nthreads) {
  t_bar->wait();

  unsigned long all_starting_sizes;  
  bool all_buckets_are_empty = false;
  typename Bucket::size_type current_bucket_size;

  while (true) {
    // first process FIX0 bucket
    current_bucket_size = buckets[FIX0_INDEX]->size();
    {
      amplusplus::scoped_epoch_value epoch(transport, current_bucket_size,
				     all_starting_sizes); 
      for (typename Bucket::size_type i = tid ; 
	   i < current_bucket_size ; i+= nthreads) {
	Vertex dest_vertex = (*buckets[FIX0_INDEX])[i];
	// inform all neighbors
	std::set<Vertex> adjacencies;
	BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
	  Vertex u = target(e, g);
	  if (logical_id(u) > logical_id(dest_vertex)) {
	    if (adjacencies.insert(u).second) {
#ifdef PRINT_DEBUG
	      std::cout << "Sending message to : " << u << " source : " << dest_vertex
			<< " mis[dest_vertex]=" << mis[dest_vertex]
			<< std::endl;
#endif
	      work_item_t wi = construct_wi(u, dest_vertex, mis[dest_vertex]);
	      relax_msg.send(wi);
	    }
	  }
	}
      }
    }

    if (tid == 0) {
      buckets[FIX0_INDEX]->clear();
    }

    t_bar->wait();

    if (all_starting_sizes == 0) {
#ifdef PRINT_DEBUG
      std::cout << "[A]All starting sizes are 0 -- breaking ...." 
		<< " FIX1 bucket size : " << buckets[FIX1_INDEX]->size()
		<< std::endl;
#endif
      break;
    } 
#ifdef PRINT_DEBUG
    else {
      std::cout << "[A]All starting sizes = " << all_starting_sizes 
		<< " FIX1 bucket size : " << buckets[FIX1_INDEX]->size()
		<< std::endl;
    }
#endif

    // done with FIX0 buckets, process FIX1 buckets
    current_bucket_size = buckets[FIX1_INDEX]->size();
    {
      amplusplus::scoped_epoch_value epoch(transport, current_bucket_size,
				     all_starting_sizes); 
      for (typename Bucket::size_type i = tid ; 
	   i < current_bucket_size ; i+= nthreads) {
	Vertex dest_vertex = (*buckets[FIX1_INDEX])[i];
#ifdef PRINT_DEBUG
	std::cout << "############# Processing state updates for " << dest_vertex << std::endl;
#endif
	std::set<Vertex> adjacencies;
	BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
	  Vertex u = target(e, g);
	  if (logical_id(u) > logical_id(dest_vertex)) { // we know that all lower neighbors are fixed, therefore we only need to send messages to higher neighbors
	    if (adjacencies.insert(u).second) {
	      work_item_t wi = construct_wi(u, dest_vertex, mis[dest_vertex]);
	      relax_msg.send(wi);
	    }
	  }
	}
      }
    }

    if (tid == 0) {
      buckets[FIX1_INDEX]->clear();
    }

    t_bar->wait();


    if (all_starting_sizes == 0) {
#ifdef PRINT_DEBUG
      std::cout << "[A]All starting sizes are 0 -- breaking ...." 
		<< " FIX0 bucket size : " << buckets[FIX0_INDEX]->size()
		<< std::endl;
#endif
      break;
    }
  }

}

template<MIS_PARMS_DELTA>
struct MIS_TYPE_DELTA::
processing_function {
  
  processing_function() : self(NULL) {}
  processing_function(maximal_independent_set_delta& self) : self(&self) {}
  
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
    self->process(data, tid);
  }

protected:
  maximal_independent_set_delta* self;
};

}}}
#endif
