// Copyright (C) 2007-2014 The Trustees of Indiana University.

#ifndef BOOST_GRAPH_KLA_SSSP_BUFFER_HPP
#define BOOST_GRAPH_KLA_SSSP_BUFFER_HPP

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
#include "boost/tuple/tuple.hpp"
//#include <tuple>
#include <boost/graph/distributed/owner_defs.hpp>

namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
         typename DistanceMap, 
         typename EdgeWeightMap, 
	 typename WorkStats,
         typename Bucket = append_buffer<std::pair<typename graph_traits<Graph>::vertex_descriptor,
						   std::pair<typename property_traits<EdgeWeightMap>::value_type, size_t > >, 10u>,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class kla_shortest_paths_buffer {
  typedef kla_shortest_paths_buffer<Graph, DistanceMap, EdgeWeightMap, Bucket> 
    self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;
  typedef typename property_traits<EdgeWeightMap>::value_type Dist;

  typedef typename std::vector<Bucket*>::size_type BucketIndex;

  typedef std::pair<Vertex, std::pair<Dist, size_t> > vertex_distance_data;
  typedef std::pair<Vertex, Dist> vertex_distance_pair;
  //typedef tuple<Vertex, Dist, size_t> vertex_distance_data_klevel;

  struct vertex_distance_handler;

  struct minimum_pair_first
  {
    template<typename T>
    const T& operator()(const T& x, const T& y) const { return x.first < y.first ? x : y; }

    template<typename F>
    struct result {
      typedef typename boost::function_traits<F>::arg1_type type;
    };
  };


  typedef typename MessageGenerator::template call_result<vertex_distance_data, vertex_distance_handler, owner_from_pair<OwnerMap, vertex_distance_data>, amplusplus::idempotent_combination_t<minimum_pair_first > >::type
    RelaxMessage;

public:

  kla_shortest_paths_buffer(Graph& g,
                            DistanceMap distance, 
                            EdgeWeightMap weight,
                            amplusplus::transport &t,                            
                            size_t k_level,
                            int offs,
			    WorkStats& stats,
                            MessageGenerator message_gen = 
                            MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_distance_data>(), 0)),
      g(g), transport(t), distance(distance), weight(weight), 
      owner(get(vertex_owner, g)), 
      level_sync(false), 
      current_level(0), 
      buckets_processed(0), 
      current_bucket(0), 
      k_level(k_level),
      core_offset(offs),
      work_stats(stats),
    relax_msg(message_gen, transport, owner_from_pair<OwnerMap, vertex_distance_data>(owner),
		amplusplus::idempotent_combination(minimum_pair_first()))
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
  void relax(const vertex_distance_data& data);

  void find_next_bucket();

  const int dummy_first_member_for_init_order; // Unused

  const Graph& g;
  amplusplus::transport& transport;
  DistanceMap distance;
  EdgeWeightMap weight;
  OwnerMap owner;
  Vertex source;
  bool level_sync;
  size_t k_level;
  int core_offset;
  WorkStats& work_stats;
  size_t k_current_level;

  // Bucket data structure. The ith bucket contains all local vertices
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

#define KLA_SHORTEST_PATHS_BUFFER_PARMS                                   \
      typename Graph, typename DistanceMap, typename EdgeWeightMap, typename WorkStats, typename Bucket, typename MessageGenerator

#define KLA_SHORTEST_PATHS_BUFFER_TYPE                                    \
      kla_shortest_paths_buffer<Graph, DistanceMap, EdgeWeightMap, WorkStats, Bucket, MessageGenerator>

template<KLA_SHORTEST_PATHS_BUFFER_PARMS>
void
KLA_SHORTEST_PATHS_BUFFER_TYPE::initialize()
{

  relax_msg.set_handler(vertex_distance_handler(*this));

  // Setup distance map
  distance.set_consistency_model(0);
  set_property_map_role(vertex_distance, distance);

  // Set the currently processing k-level
  k_current_level = k_level;

  // Initialize buckets data structure
  //

  // Extra bucket is so we don't try to insert into the bucket we're processing
  // when the index wraps
  num_buckets = 2;

  // Declare bucket data structure and index variable
  buckets.resize(num_buckets);
  for (BucketIndex i = 0 ; i < buckets.size() ; ++i) {
    shared_ptr<Bucket> p(new Bucket);
    buckets[i].swap(p);
  }

  // Initialize distance labels
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    put(distance, v, (std::numeric_limits<Dist>::max)()); 
  }
}

template<KLA_SHORTEST_PATHS_BUFFER_PARMS>
void
KLA_SHORTEST_PATHS_BUFFER_TYPE::run(Vertex s, int tid)
{
  //int debugwait = 1;
  //if(transport.rank()==0) 
  //if(tid==0)
  //while (debugwait) ;
  
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

    start_time = get_time();
    // Push the source onto the bucket
    { 
      amplusplus::scoped_epoch epoch(transport); 

      if (get(owner, s) == transport.rank() && tid == 0) {
	//	buckets[current_bucket]->push_back(vertex_distance_data(s,std::make_pair(0,0)));
	relax(vertex_distance_data(s,std::make_pair(0,0)));
      }
    }

    t_bar->wait();
    
    typename Bucket::size_type current_bucket_start, current_bucket_end;
    
    int count_epoch = 0;
    do {
      current_bucket_start = 0;
      current_bucket_end = buckets[current_bucket]->size();
      // Process current bucket
      do {
	unsigned long all_starting_sizes;
	const unsigned long starting_size = current_bucket_end - current_bucket_start;
        
	t_bar->wait();
	{
	  amplusplus::scoped_epoch_value epoch(transport, starting_size, all_starting_sizes);
#ifdef PRINT_DEBUG
	  if (tid == 0)
	    std::cerr << transport.rank() << ": processing edges in bucket "
		      << current_bucket << " [" << current_bucket_start << ", " << current_bucket_end << ")\n";
#endif
	  
	  while (current_bucket_start != current_bucket_end) {
	    for (typename Bucket::size_type i = current_bucket_start + tid ; i < current_bucket_end ; i+= nthreads) {
	      Vertex v = (*buckets[current_bucket])[i].first;
	      Dist dist = (*buckets[current_bucket])[i].second.first;
	      size_t v_level = (*buckets[current_bucket])[i].second.second;
	      Dist dv = get(distance, v);
	    
	      //std::cout<<"\nProcessing vertex:  "<<v <<" with distance "<<dv<<" and checking whether we got a better distance "<<dist<<std::endl;
	      //if we already have a better distance in the distance map
              if (dv < dist) {
		continue;
              }
	      //Otherwise we got  a better distance.  Relax edges
	      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
		Vertex u = target(e, g);
		Dist we = get(weight, e);
		vertex_distance_data new_vd(u, std::make_pair(dist + we, v_level+1));
		relax_msg.send(new_vd);
	      }
	    }
	  
	    // Wait for all threads to finish current bucket
	    t_bar->wait();
	    current_bucket_start = current_bucket_end;
	    current_bucket_end = buckets[current_bucket]->size();       
	    t_bar->wait();

	    //if (level_sync) break;
	  }
	}
	count_epoch += 1;
	// If all processes are done with the current bucket:
	// 1. clear the current bucket
	// 2. increase current level
	// 3. find the next bucket to work on  
	if (all_starting_sizes == 0) {
	  // current bucket is done in all processes
	  if (tid == 0) buckets[current_bucket]->clear();
	  t_bar->wait();
	  
#ifdef PRINT_DEBUG
	  std::cerr << tid << "@" << transport.rank() << ": Current bucket " << current_bucket
		    << " size " << buckets[current_bucket]->size() << std::endl;
#endif
	  assert(buckets[current_bucket]->size() == 0);
         
	  // find next bucket with work and update current_level 
	  BucketIndex old_bucket = current_bucket;
          
	  t_bar->wait();
	  if (tid == 0) {
	    find_next_bucket();

	    if (current_bucket != (std::numeric_limits<BucketIndex>::max)()) {
	      current_level += current_bucket - old_bucket;
	      ++buckets_processed;
	    }
	    
	    // update k_current_level
	    k_current_level += k_level;
	  }
	  t_bar->wait();
      
	  break; // exit the do..while loop for the current bucket
	} else // Update the end of the bucket and continue processing it
	  current_bucket_end = buckets[current_bucket]->size();
      }while(true);
     
      // If there are no non-empty buckets in any process, we're done
    } while(current_bucket != (std::numeric_limits<BucketIndex>::max)());
    
  }
}

template<KLA_SHORTEST_PATHS_BUFFER_PARMS>
void
KLA_SHORTEST_PATHS_BUFFER_TYPE::relax(const vertex_distance_data& data) 
{

#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_edges(amplusplus::detail::get_thread_id());
#endif  

  Vertex v = data.first;
  Dist d = data.second.first;
  size_t k_vertex = data.second.second;
  using boost::parallel::val_compare_and_swap;

  Dist old_dist = get(distance, v), last_old_dist;
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

      // We got a better distance. If current
      // k-level is less than or equal to vertex k-level
      // straight a way relax it
      if (k_vertex <= k_current_level) {
	assert((k_current_level-k_level) < k_vertex <= k_current_level);

	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  Dist we = get(weight, e);
	  vertex_distance_data new_vd(u, std::make_pair(d + we, k_vertex+1));
	  relax_msg.send(new_vd);
	}
      } else { // otherwise we put vertex to next bucket
	buckets[(current_bucket+1)%2]->push_back(data);
      }

      break; // No need to insert into current bucket now
    }
  }
#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(amplusplus::detail::get_thread_id());
#endif
  return;
}

template<KLA_SHORTEST_PATHS_BUFFER_PARMS>
void
KLA_SHORTEST_PATHS_BUFFER_TYPE::find_next_bucket()
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

template<KLA_SHORTEST_PATHS_BUFFER_PARMS>
struct KLA_SHORTEST_PATHS_BUFFER_TYPE::
vertex_distance_handler {
  
  vertex_distance_handler() : self(NULL) {}
  vertex_distance_handler(kla_shortest_paths_buffer& self) : self(&self) {}
  
  void operator() (const vertex_distance_data& data) const {
    self->relax(data);
  }

protected:
  kla_shortest_paths_buffer* self;
};


} } } // end namespace boost::graph::distributed

#endif // BOOST_GRAPH_KLA_SHORTEST_PATHS_BUFFER_HPP
