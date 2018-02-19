// Copyright (C) 2007-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

/**************************************************************************
 * This source file implements the Delta-stepping algorithm:              *
 *                                                                        *
 *   Ulrich Meyer and Peter Sanders. Parallel Shortest Path for Arbitrary *
 *   Graphs. In Proceedings from the 6th International Euro-Par           *
 *   Conference on Parallel Processing, pages 461--470, 2000.             *
 *                                                                        * 
 *   Ulrich Meyer, Peter Sanders: [Delta]-stepping: A Parallelizable      *
 *   Shortest Path Algorithm. J. Algorithms 49(1): 114-152, 2003.         *
 *                                                                        *
 * There are several potential optimizations we could still perform for   *
 * this algorithm implementation:                                         *
 *                                                                        *
 *   - Implement "shortcuts", which bound the number of reinsertions      *
 *     in a single bucket (to one). The computation of shortcuts looks    *
 *     expensive in a distributed-memory setting, but it could be         *
 *     ammortized over many queries.                                      *
 *                                                                        *
 *   - The size of the "buckets" data structure can be bounded to         *
 *     max_edge_weight/delta buckets, if we treat the buckets as a        *
 *     circular array.                                                    *
 *                                                                        *
 *   - If we partition light/heavy edges ahead of time, we could improve  *
 *     relaxation performance by only iterating over the right portion    *
 *     of the out-edge list each time.                                    *
 *                                                                        *
 *   - Implement a more sophisticated algorithm for guessing "delta",     *
 *     based on the shortcut-finding algorithm.                           *
 **************************************************************************/
#ifndef BOOST_GRAPH_DELTA_STEPPING_SHORTEST_PATHS_HPP
#define BOOST_GRAPH_DELTA_STEPPING_SHORTEST_PATHS_HPP

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
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <algorithm> // for std::min, std::max
#include <boost/format.hpp>
#include <iostream>
#include <atomic>

template <typename OwnerMap, typename VertexDistanceData>
struct vertex_distance_owner {
  explicit vertex_distance_owner(const OwnerMap& owner) : owner(owner) {}

  const OwnerMap& owner;
};

template <typename OwnerMap, typename VertexDistanceData>
typename boost::property_traits<OwnerMap>::value_type
get(const vertex_distance_owner<OwnerMap, VertexDistanceData>& o, const VertexDistanceData& data)
{ return get(o.owner, data.first); }

namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
	 typename DistanceMap, 
	 typename EdgeWeightMap, 
	 typename WorkStats,
         typename Bucket = append_buffer<typename graph_traits<Graph>::vertex_descriptor, 10u>,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class delta_stepping_shortest_paths {
  typedef delta_stepping_shortest_paths<Graph, DistanceMap, EdgeWeightMap, Bucket> 
    self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;
  typedef typename property_traits<EdgeWeightMap>::value_type Dist;

  typedef typename std::vector<Bucket*>::size_type BucketIndex;

  typedef std::pair<Vertex, Dist> vertex_distance_data;

  struct vertex_distance_handler;

  typedef typename MessageGenerator::template call_result<vertex_distance_data, vertex_distance_handler, 
				     vertex_distance_owner<OwnerMap, vertex_distance_data>,
				     amplusplus::idempotent_combination_t<boost::parallel::minimum<Dist>,
									  Dist> >::type
    RelaxMessage;


public:

  delta_stepping_shortest_paths(Graph& g,
                                DistanceMap& distance, 
                                EdgeWeightMap& weight,
				amplusplus::transport &t,
                                Dist delta,
				int offs,
				WorkStats& stats,
                                MessageGenerator message_gen = 
                                  MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_distance_data>(), 0)),
      g(g), transport(t), distance(distance), weight(weight), delta(delta), 
    core_offset(offs),
    work_stats(stats),
    owner(get(vertex_owner, g)), 
      level_sync(false), current_level(0), buckets_processed(0), current_bucket(0),
      relax_msg(message_gen, transport, vertex_distance_owner<OwnerMap, vertex_distance_data>(owner),
		amplusplus::idempotent_combination(boost::parallel::minimum<Dist>(), std::numeric_limits<Dist>::max()))
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
  void relax(Vertex v, Dist x);

  void find_next_bucket();

  const int dummy_first_member_for_init_order; // Unused

  const Graph& g;
  amplusplus::transport& transport;
  DistanceMap& distance;
  EdgeWeightMap& weight;
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
};

#define DELTA_STEPPING_SHORTEST_PATHS_PARMS                                   \
      typename Graph, typename DistanceMap, typename EdgeWeightMap, typename WorkStats, typename Bucket, typename MessageGenerator

#define DELTA_STEPPING_SHORTEST_PATHS_TYPE                                    \
      delta_stepping_shortest_paths<Graph, DistanceMap, EdgeWeightMap, WorkStats, Bucket, MessageGenerator>

template<DELTA_STEPPING_SHORTEST_PATHS_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_TYPE::initialize()
{

  relax_msg.set_handler(vertex_distance_handler(*this));

  // Setup distance map
  distance.set_consistency_model(0);
  set_property_map_role(vertex_distance, distance);

  using boost::parallel::all_reduce;
  using boost::parallel::maximum;
  using std::max;

  // Compute the maximum edge weight
  Dist max_edge_weight = 0;
  Dist local_max = 0;
  // CULPRIT !
  BGL_FORALL_VERTICES_T(u, g, Graph) {
    BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
      auto w = get(weight, e);
      if (local_max < w)
	local_max = w;
    }
  }

  // do an all reduce to find the maximum
  MPI_Allreduce(&local_max, &max_edge_weight, 1, 
		MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);

  //  all_reduce<Dist, maximum<Dist> > r(transport, maximum<Dist>());
  //max_edge_weight = r(max_edge_weight);

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

  // Declare bucket data structure and index variable
  buckets.resize(num_buckets);
  for (BucketIndex i = 0 ; i < buckets.size() ; ++i) {
    shared_ptr<Bucket> p(new Bucket);
    buckets[i].swap(p);
  }

  shared_ptr<Bucket> p(new Bucket);
  deleted_vertices.swap(p);

  // Initialize distance labels
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    put(distance, v, (std::numeric_limits<Dist>::max)()); 
  }
}

template<DELTA_STEPPING_SHORTEST_PATHS_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_TYPE::run(Vertex s, int tid)
{
  /**
   *  NOTE: We never remove a vertex from a higher bucket when
   *  placing it in a lower one because thread-safe insertion is
   *  cheap as long as we don't have to handle removal
   *  simultaneously, this will result in additional work and higher
   *  memory consumption, but the alternative is locking which is
   *  sure to be slow.
   */
  int count_epoch = 0 ; 
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    int nthreads = transport.get_nthreads();
    
    if (tid == 0)
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    
    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized 
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

    // last passing thread will update the correct start time
    start_time = get_time();
    // Push the source onto the queue
    { 
      amplusplus::scoped_epoch epoch(transport); 

      if (get(owner, s) == transport.rank() && tid == 0)
        relax(s, 0);
    }    

    t_bar->wait();
    
    typename Bucket::size_type current_bucket_start, current_bucket_end;
    
    do {
      
      current_bucket_start = 0;
      current_bucket_end = buckets[current_bucket]->size();
      
      // Process current bucket
      do {
        
        unsigned long all_starting_sizes;
        const unsigned long starting_size = current_bucket_end - current_bucket_start;
        count_epoch = 0;
        t_bar->wait();
        {
          amplusplus::scoped_epoch_value epoch(transport, starting_size, all_starting_sizes);
          
#ifdef PRINT_DEBUG
          if (tid == 0)
            std::cerr << transport.rank() << ": processing light edges in bucket "
                      << current_bucket << " [" << current_bucket_start << ", " << current_bucket_end << ")\n";
#endif
	  
	  while (current_bucket_start != current_bucket_end) {
	    
	    for (typename Bucket::size_type i = current_bucket_start + tid ; i < current_bucket_end ; i+= nthreads) {
	      Vertex v = (*buckets[current_bucket])[i];
	      Dist dv = distance[v];
	      
              // This is an optimization that is not specified in the algorithm
	      if (dv < delta * (Dist)current_level) {
		continue;
              }

	      // Add v to set of deleted vertices
	      deleted_vertices->push_back(v);
	      
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
              }
            }
            
            // Wait for all threads to finish current bucket
            t_bar->wait();
            current_bucket_start = current_bucket_end;
            current_bucket_end = buckets[current_bucket]->size();       
            t_bar->wait();

	    if (level_sync) break;
	  } // No more re-insertions into this bucket to process
	}
	count_epoch += 1;
	// If all processes are done with the current bucket:
	//   1. clear the current bucket
	//   2. process heavy edges
	//   3. find the next bucket to work on 
	if (all_starting_sizes == 0) {
	  
	  if (tid == 0) buckets[current_bucket]->clear();
	  t_bar->wait();
	  
	  // Process heavy edges
	  // TODO: If we had a separate transport here we could do single-level termination detection
	  {
	    amplusplus::scoped_epoch epoch(transport);
	    
	    int nthreads = transport.get_nthreads();
	    for (typename Bucket::size_type i = tid ; i < deleted_vertices->size() ; i += nthreads) {
	      Vertex v = (*deleted_vertices)[i];
	      Dist dv = distance[v];
	      
	      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
		Vertex u = target(e, g);
		Dist we = get(weight, e);
		if (we > delta) // heavy edge
		  relax_msg.send(vertex_distance_data(u, dv + we));
	      }
	    }
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
  
          break; // exit the do..while loop for the current bucket
        } else // Update the end of the bucket and continue processing it
          current_bucket_end = buckets[current_bucket]->size();
        
      } while(true);
      
      // If there are no non-empty buckets in any process, we're done
    } while(current_bucket != (std::numeric_limits<BucketIndex>::max)());
  }
}

template<DELTA_STEPPING_SHORTEST_PATHS_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_TYPE::relax(Vertex v, Dist d) 
{
  using boost::parallel::val_compare_and_swap;

  Dist old_dist = distance[v], last_old_dist;
  while (d < old_dist) {
    last_old_dist = old_dist;
    old_dist = val_compare_and_swap(&distance[v], old_dist, d);
    if (last_old_dist == old_dist) {
#ifdef PBGL2_PRINT_WORK_STATS
      int tid = amplusplus::detail::get_thread_id();
      if(old_dist < std::numeric_limits<Dist>::max()) { 
	work_stats.increment_invalidated(tid);
      } else {
	work_stats.increment_useful(tid);
      }
#endif
      // Insert vertex into new bucket, note we don't remove it from any other
      // buckets it might be in 
      BucketIndex new_index = 
        static_cast<BucketIndex>((d - (current_level * delta)) / delta);

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
	
	deleted_vertices->push_back(v);

        break; // No need to insert into current bucket now
      }

      new_index = (current_bucket + new_index) % num_buckets;
      buckets[new_index]->push_back(v);

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

template<DELTA_STEPPING_SHORTEST_PATHS_PARMS>
void
DELTA_STEPPING_SHORTEST_PATHS_TYPE::find_next_bucket()
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

template<DELTA_STEPPING_SHORTEST_PATHS_PARMS>
struct DELTA_STEPPING_SHORTEST_PATHS_TYPE::
vertex_distance_handler {
  
  vertex_distance_handler() : self(NULL) {}
  vertex_distance_handler(delta_stepping_shortest_paths& self) : self(&self) {}
  
  void operator() (const vertex_distance_data& data) const {
    self->work_stats.increment_edges(amplusplus::detail::get_thread_id());

    if (data.second < get(self->distance, data.first))
      self->relax(data.first, data.second);
#ifdef PBGL2_PRINT_WORK_STATS
    else
      self->work_stats.increment_rejected(amplusplus::detail::get_thread_id());
#endif
  }

protected:
  delta_stepping_shortest_paths* self;
};


} } } // end namespace boost::graph::distributed

#endif // BOOST_GRAPH_DELTA_STEPPING_SHORTEST_PATHS_HPP
