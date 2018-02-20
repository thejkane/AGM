// Copyright (C) 2007-2018 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_DELTA_STEPPING_SHORTEST_PATHS_NODE_HPP
#define BOOST_GRAPH_DELTA_STEPPING_SHORTEST_PATHS_NODE_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/graph/distributed/priority_q_defs.hpp>
#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <algorithm> // for std::min, std::max
#include <boost/format.hpp>
#include <iostream>
#include <atomic>

namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
         typename DistanceMap, 
         typename EdgeWeightMap, 
	 typename WorkStats,
	 typename PriorityQueueGenerator = node_priority_queue_gen,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class delta_stepping_shortest_paths_node {
  typedef delta_stepping_shortest_paths_node<Graph, DistanceMap, EdgeWeightMap, PriorityQueueGenerator> 
    self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;
  typedef typename property_traits<EdgeWeightMap>::value_type Dist;
  typedef std::pair<Vertex, Dist> vertex_distance_data;

  // Constructing bucket type ...
  struct default_comparer {
    bool operator()(const vertex_distance_data& vd1, const vertex_distance_data& vd2) const {
      return vd1.second > vd2.second;
    }
  };

  typedef typename PriorityQueueGenerator::template queue<vertex_distance_data, 
							  default_comparer>::type Bucket;

  // Bucket is an ordering container
  //  typedef typename PriorityQueueType Bucket;

  typedef typename std::vector<Bucket*>::size_type BucketIndex;

  struct vertex_distance_handler;

  typedef typename MessageGenerator::template call_result<vertex_distance_data, vertex_distance_handler, 
				     vertex_distance_owner<OwnerMap, vertex_distance_data>,
				     amplusplus::idempotent_combination_t<boost::parallel::minimum<Dist>,
									  Dist> >::type
    RelaxMessage;


public:

  delta_stepping_shortest_paths_node(Graph& g,
				     DistanceMap distance, 
				     EdgeWeightMap weight,
				     amplusplus::transport &t,
				     Dist delta,
				     int offs,
                                     WorkStats& stats,
				     int freq,
				     MessageGenerator message_gen = 
				     MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_distance_data>(), 0)),
      g(g), transport(t), distance(distance), weight(weight), delta(delta), 
      core_offset(offs),
      work_stats(stats),
      owner(get(vertex_owner, g)), 
      level_sync(false), current_level(0), buckets_processed(0), current_bucket(0),
      relax_msg(message_gen, transport, vertex_distance_owner<OwnerMap, vertex_distance_data>(owner),
		amplusplus::idempotent_combination(boost::parallel::minimum<Dist>(), std::numeric_limits<Dist>::max())),
            flushFrequency(freq)
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

  void set_source(Vertex s) { source = s; } // for threaded execution
  void set_level_sync() { level_sync = true; } // force level-synchronized exploration
  void operator() (int tid) { run(source, tid); }

  void run(Vertex s, int tid = 0);

  BucketIndex get_num_levels() { return buckets_processed; }
  time_type get_start_time() { return start_time; }

protected:
  void initialize();
  
  // Relax the edge (u, v), creating a new best path of distance x.
  void relax(const vertex_distance_data& vd);

  void find_next_bucket();

  const int dummy_first_member_for_init_order; // Unused

  const Graph& g;
  amplusplus::transport& transport;
  DistanceMap distance;
  EdgeWeightMap weight;
  Dist delta;
  int core_offset;
  WorkStats& work_stats;
  OwnerMap owner;
  Vertex source;
  bool level_sync;

  // Bucket data structure. The ith bucket contains all local vertices
  // with (tentative) distance in the range [i*delta,
  // (i+1)*delta). 
  std::vector<shared_ptr<Bucket> > buckets;

  // Bucket to hold vertices deleted at each level
  shared_ptr<Bucket> deleted_vertices;
  BucketIndex current_level; // How many buckets have we processed?
  BucketIndex buckets_processed; // Stats tracking
  BucketIndex num_buckets;

  // Shared thread state to make sure we're all on the same page
  BucketIndex current_bucket;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  RelaxMessage relax_msg;
  time_type start_time;
  int flushFrequency;
};

#define DELTA_STEPPING_SHORTEST_PATHS_NODE_PARMS                                   \
      typename Graph, typename DistanceMap, typename EdgeWeightMap, typename WorkStats, typename PriorityQueueGenerator, typename MessageGenerator

#define DELTA_STEPPING_SHORTEST_PATHS_NODE_TYPE                                    \
      delta_stepping_shortest_paths_node<Graph, DistanceMap, EdgeWeightMap, WorkStats, PriorityQueueGenerator, MessageGenerator>

template<DELTA_STEPPING_SHORTEST_PATHS_NODE_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_NODE_TYPE::initialize()
{
  int nthreads = transport.get_nthreads();

  relax_msg.set_handler(vertex_distance_handler(*this));

  // Setup distance map
  distance.set_consistency_model(0);
  set_property_map_role(vertex_distance, distance);

  using boost::parallel::all_reduce;
  using boost::parallel::maximum;
  using std::max;

  // Compute the maximum edge weight
  Dist max_edge_weight = 0;

  BGL_FORALL_VERTICES_T(u, g, Graph) {
    BGL_FORALL_OUTEDGES_T(u, e, g, Graph)
      max_edge_weight = max BOOST_PREVENT_MACRO_SUBSTITUTION (max_edge_weight, get(weight, e));
  }

  all_reduce<Dist, maximum<Dist> > r(transport, maximum<Dist>());
  max_edge_weight = r(max_edge_weight);

  // If delta wasn't supplied in the ctor initialize it
  if (delta == 0) {

    // Compute the maximum edge degree
    Degree max_degree = 0;
    BGL_FORALL_VERTICES_T(u, g, Graph) {
      max_degree = max BOOST_PREVENT_MACRO_SUBSTITUTION (max_degree, out_degree(u, g));
    }
    // max_degree = all_reduce(process_group(g), max_degree, maximum<Degree>());
    all_reduce<Degree, maximum<Degree> > r(transport, maximum<Degree>());
    max_degree = r(max_degree);
    
    // Take a guess at delta, based on what works well for random
    // graphs.
    delta = max_edge_weight / max_degree;
    if (delta == 0)
      delta = 1;
  }

  //
  // Initialize buckets data structure
  //

  // Extra bucket is so we don't try to insert into the bucket we're processing
  // when the index wraps
  num_buckets = max_edge_weight / delta;
  num_buckets += (num_buckets * delta < (BucketIndex)max_edge_weight) ? 2 : 1;

  // + 1 is for deleted_vertices bucket
  unsigned long capacity_per_bucket = (DEFAULT_HEAP_SIZE / (num_buckets + 1));

  // Declare bucket data structure and index variable
  buckets.resize(num_buckets);
  for (BucketIndex i = 0 ; i < buckets.size() ; ++i) {
    shared_ptr<Bucket> p(new Bucket(nthreads, capacity_per_bucket));
    buckets[i].swap(p);
  }

  shared_ptr<Bucket> p(new Bucket(nthreads, capacity_per_bucket));
  deleted_vertices.swap(p);

  // Initialize distance labels
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    put(distance, v, (std::numeric_limits<Dist>::max)()); 
  }
}

template<DELTA_STEPPING_SHORTEST_PATHS_NODE_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_NODE_TYPE::run(Vertex s, int tid) {
  /**
   *  NOTE: We never remove a vertex from a higher bucket when
   *  placing it in a lower one because thread-safe insertion is
   *  cheap as long as we don't have to handle removal
   *  simultaneously, this will result in additional work and higher
   *  memory consumption, but the alternative is locking which is
   *  sure to be slow.
   */
  int doFlushCounter = 0;
  int count_epoch = 0 ; 
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    int nthreads = transport.get_nthreads();

    if (tid == 0)
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    
    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized 
    { amplusplus::scoped_epoch epoch(transport); } 

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

#ifdef LIB_CDS
    // Attach thread to CDS; main thread is already attached
    if (tid != 0)
      cds::threading::Manager::attachThread();
#endif

    start_time = get_time();
    // Push the source onto the queue
    { 
      amplusplus::scoped_epoch epoch(transport); 

      if (get(owner, s) == transport.rank() && tid == 0)
        relax(vertex_distance_data(s, 0));
    }    

    t_bar->wait();

    while(current_bucket != (std::numeric_limits<BucketIndex>::max)()) { // do for all buckets

      //std::cout << transport.rank() 
      //	<< " processing bucket " << current_bucket << std::endl;

      unsigned long all_process_bucket_empty = 1;
      unsigned long p_bucket_full = 0; 
      count_epoch = 0;
      // process current bucket
      while(all_process_bucket_empty != 0) {

	// wait till all threads reach here
	// before checking whether queue is empty we need to make sure
	// none of the threads are pushing elements to the queue.
	// Therefore make sure all threads reach here before checking whether queue is
	// empty
	t_bar->wait();
	if (tid == 0) { // this is only valid if we are using a node specific data structure
	  if (!buckets[current_bucket]->empty())
	    p_bucket_full = 1;
	  else
	    p_bucket_full = 0;
	}
	//t_bar->wait(); // TODO Not sure whether we really need this

	{
	  //	  std::cout << "p_bucket_not_empty : " << p_bucket_not_empty << std::endl;
	  amplusplus::scoped_epoch_value epoch(transport, 
					       p_bucket_full, 
					       all_process_bucket_empty);

 	  //std::cout << transport.rank() 
	  //    << " all_process_bucket_empty : " 
	  //    << all_process_bucket_empty << std::endl;

	  vertex_distance_data vd;
	  while(buckets[current_bucket]->pop(vd, tid)) { // extract an item from the current bucket
	    
	    assert(all_process_bucket_empty != 0);
	    // process light edges
#ifdef PRINT_DEBUG
	    if (tid == 0)
	      std::cerr << transport.rank() << ": processing light edges in bucket "
			<< current_bucket << std::endl;
#endif
	    Vertex v = vd.first;
	    Dist dv = get(distance, v);

	    // This is an optimization that is not specified in the algorithm
	    if (dv < delta * (Dist)current_level) {
	      continue;
            }

	    // Add v to set of deleted vertices
	    deleted_vertices->put(vd, tid);
	      
	    // Relax light edges
	    BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	      Vertex u = target(e, g);
	      Dist we = get(weight, e);
	      if (we <= delta) {// light edge
		relax_msg.send(vertex_distance_data(u, dv + we));
#ifdef PRINT_DEBUG
		std::cerr << "  Relax " << get(get(vertex_local, g), v) << "@" 
			  << get(get(vertex_owner, g), v) << "->"
			  << get(get(vertex_local, g), u) << "@" 
			  << get(get(vertex_owner, g), u) << "  weight = "
			  << we << std::endl;
#endif
	      }
	    } // end of BGL for

	    doFlushCounter++;
	    if(doFlushCounter == flushFrequency) {
	      doFlushCounter = 0;
	      transport.get_scheduler().run_one();
	    }
	  } // end of while (bucket)
	  // TODO need to add if (level_sync) break; code; not sure still valid
	} // end of scoped epoch      
	count_epoch += 1;
      } // all processors should be done with current bucket

      assert(buckets[current_bucket]->empty());

      // If all processes are done with the current bucket:
      //   1. clear the current bucket
      //   2. process heavy edges
      //   3. find the next bucket to work on 
      if (tid == 0) buckets[current_bucket]->clear();
      t_bar->wait();
	  
      // Process heavy edges
      // TODO: If we had a separate transport here we could do single-level termination detection
      {
	amplusplus::scoped_epoch epoch(transport);
	    
	vertex_distance_data deleted_vd;
	while(deleted_vertices->pop(deleted_vd, tid)) {
	  Vertex v = deleted_vd.first;
	  Dist dv = get(distance, v);
	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    Vertex u = target(e, g);
	    Dist we = get(weight, e);
	    if (we > delta) // heavy edge
	      relax_msg.send(vertex_distance_data(u, dv + we));
	  } // end BGL for
	  
	  doFlushCounter++;
	  if(doFlushCounter == flushFrequency) {
	    doFlushCounter = 0;
	    transport.get_scheduler().run_one();
	  }
	} // end while
      } // end_epoch does a thread barrier
	  
#ifdef PRINT_DEBUG
      std::cerr << tid << "@" << transport.rank() << ": Current bucket " << current_bucket
		<< " size " << buckets[current_bucket]->size() << std::endl;
#endif
      assert(buckets[current_bucket]->size() == 0);

      if (tid == 0) deleted_vertices->clear();
          
      // find next bucket with work and update current_level 
      BucketIndex old_bucket = current_bucket;
          
      t_bar->wait();
      if (tid == 0) {
	find_next_bucket();

	if (current_bucket != (std::numeric_limits<BucketIndex>::max)()) {
	  current_level += current_bucket - old_bucket;
	  ++buckets_processed;
	}
      }
      t_bar->wait();
    } // end of while global bucket check 

#ifdef LIB_CDS
    if (tid != 0)
      cds::threading::Manager::detachThread();
#endif 

  } // end of AMPP THREAD ID
} //end of run

template<DELTA_STEPPING_SHORTEST_PATHS_NODE_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_NODE_TYPE::relax(const vertex_distance_data& vd) 
{
  Vertex v = vd.first;
  Dist d = vd.second;
  using boost::parallel::val_compare_and_swap;

  Dist old_dist = get(distance, v), last_old_dist;
  while (d < old_dist) {
    last_old_dist = old_dist;
    old_dist = val_compare_and_swap(&distance[v], old_dist, d);
    if (last_old_dist == old_dist) {
#ifdef PBGL2_PRINT_WORK_STATS
      if(old_dist < std::numeric_limits<Dist>::max()) 
	work_stats.increment_invalidated(amplusplus::detail::get_thread_id());
      else {
	work_stats.increment_useful(amplusplus::detail::get_thread_id());	
      }
#endif
      // Insert vertex into new bucket, note we don't remove it from any other
      // buckets it might be in 
      BucketIndex new_index = 
        static_cast<BucketIndex>((d - (current_level * delta)) / delta);

      int tid = amplusplus::detail::get_thread_id();
      //
      // If new_index == 0, relax all light edges and add it to deleted_vertices
      //
      if (new_index == 0) { // This vertex is bound for the current bucket 
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  Dist we = get(weight, e);
	  if (we <= delta) // light edge
	    relax_msg.send(vertex_distance_data(u, d + we));
	}
	
	deleted_vertices->put(vd, tid);

        break; // No need to insert into current bucket now
      }

      new_index = (current_bucket + new_index) % num_buckets;
      buckets[new_index]->put(vd, tid);

      // new_distance now == get(distance, v) so we exit automatically
      // but this saves us the conditional test
      break;
    }
  }
#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(amplusplus::detail::get_thread_id());
#endif
  return;
}

template<DELTA_STEPPING_SHORTEST_PATHS_NODE_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_NODE_TYPE::find_next_bucket()
{
  using boost::parallel::all_reduce;
  using boost::parallel::minimum;

  BucketIndex old_bucket = current_bucket;
  BucketIndex max_bucket = (std::numeric_limits<BucketIndex>::max)();

  current_bucket = (current_bucket + 1) % buckets.size();
  while (current_bucket != old_bucket && buckets[current_bucket]->empty())
    current_bucket = (current_bucket + 1) % buckets.size();
  
  if (current_bucket == old_bucket) 
    current_bucket = max_bucket;
  
  // If we wrapped, project index past end of buckets to use min()
  if (current_bucket < old_bucket) current_bucket += buckets.size();
  
  all_reduce<BucketIndex, minimum<BucketIndex> > r(transport, minimum<BucketIndex>());
  current_bucket = r(current_bucket);
  
  // Map index back into range of buckets
  if (current_bucket != max_bucket) 
    current_bucket %= buckets.size();
}

template<DELTA_STEPPING_SHORTEST_PATHS_NODE_PARMS>
struct DELTA_STEPPING_SHORTEST_PATHS_NODE_TYPE::
vertex_distance_handler {
  
  vertex_distance_handler() : self(NULL) {}
  vertex_distance_handler(delta_stepping_shortest_paths_node& self) : self(&self) {}
  
  void operator() (const vertex_distance_data& data) const
  {
#ifdef PBGL2_PRINT_WORK_STATS
    self->work_stats.increment_edges(amplusplus::detail::get_thread_id());
#endif

    if (data.second < get(self->distance, data.first))
      self->relax(data);
#ifdef PBGL2_PRINT_WORK_STATS
    else
      self->work_stats.increment_rejected(amplusplus::detail::get_thread_id());
#endif
  }

protected:
  delta_stepping_shortest_paths_node* self;
};


} } } // end namespace boost::graph::distributed

#endif // BOOST_GRAPH_DELTA_STEPPING_SHORTEST_PATHS_HPP
