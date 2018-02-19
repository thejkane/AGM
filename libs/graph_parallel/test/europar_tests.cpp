// Copyright (C) 2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine

// A small test of nearly everything

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <inttypes.h>
#include <cstdlib>

//#define DISABLE_SELF_SEND_CHECK
#define AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
#define PBGL2_PRINT_WORK_STATS
// #define PRINT_STATS
// #define PRINT_DEBUG
#define AMPLUSPLUS_PRINT_HIT_RATES
#define BFS_SV_CC_HACK
// #define DS_SHARED_REDUCTIONS
// #define BFS_VECTOR
#define DISABLE_BFS_BITMAP
#ifdef __bgp__
#define BGP_REPORT_MEMORY 1 // if >1 all ranks will report memory usage
#endif
#define NUMBER_COMPONENTS
// #define PBGL_TIMING // turn on local component counting, sequential timing,  and such in CC
// #define PRINT_ET
#define BARRIER_TRANS
// #define CLONE

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>
#include <am++/counter_coalesced_message_type.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>

#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/recursive_rmat_generator.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/graph500_generator.hpp>
#include <boost/graph/permute_graph.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <boost/graph/distributed/self_stabilizing_shortest_paths.hpp>
#include <boost/graph/distributed/connected_components.hpp>
#include <boost/graph/distributed/page_rank.hpp>
#include <boost/graph/relax.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>

#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/distributed/distributed_control.hpp>
#define HAS_NOT_BEEN_SEEN

// Testing purpose
#define DISABLE_SELF_SEND_CHECK


int _RANK;

using namespace boost;

//
// Memory utilization on BG/P
//
#ifdef BGP_REPORT_MEMORY
/* compile with -I/bgsys/drivers/ppcfloor/arch/include */

#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <malloc.h>

/* gives total memory size per MPI process */
void
getMemSize(long long *mem)
{
    long long total;
    int node_config;
    _BGP_Personality_t personality;

    Kernel_GetPersonality(&personality, sizeof(personality));
    total = BGP_Personality_DDRSizeMB(&personality);

    node_config  = BGP_Personality_processConfig(&personality);
    if (node_config == _BGP_PERS_PROCESSCONFIG_VNM) total /= 4;
    else if (node_config == _BGP_PERS_PROCESSCONFIG_2x2) total /= 2;
    total *= 1024*1024;

    *mem = total;
}

/* gives total memory used 
 *  NOTE: this does not account for static data, only heap
 */
void
getUsedMem(long long *mem)
{
    long long alloc;
    struct mallinfo m;

    m = mallinfo();
    alloc = m.hblkhd + m.uordblks;

    *mem = alloc;
}

/* gives available memory 
 *  NOTE: this does not account for static data
 */
void
getFreeMem(long long *mem)
{
    long long total, alloc;

    getMemSize(&total);
    getUsedMem(&alloc);

    *mem = total - alloc;
}
#endif // BGP_REPORT_MEMORY

//
// Timing code
//

#include <time.h>
#include <sys/time.h>

typedef double time_type;

inline time_type get_time()
{
  return MPI_Wtime();
#if 0
  timeval tp;
  gettimeofday(&tp, 0);
  return tp.tv_sec + tp.tv_usec / 1000000.0;
#endif
}

std::string print_time(time_type t)
{
  std::ostringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
  return out.str();
}

//
// Performance counters
//
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS

std::vector<time_type> epoch_times;
// These should really be atomic
typedef unsigned long long flush_type;
typedef amplusplus::detail::atomic<flush_type> atomic_flush_type;
#define FLUSH_MPI_TYPE MPI_UNSIGNED_LONG
atomic_flush_type *flushes;
flush_type flushes_size;
std::vector<flush_type> all_flushes, cumulative_flushes;
atomic_flush_type full{}, messages_received{};
flush_type all_full{}, cumulative_full{}, all_messages{}, cumulative_messages{};

void print_and_clear_epoch_times()
{
  if(_RANK == 0) {
    std::cout << "There were " << epoch_times.size() << " epochs." << std::endl;
    std::cout << "Epoch times ";
    BOOST_FOREACH(time_type t, epoch_times) {
      std::cout << print_time(t) << " ";
    }
    std::cout << "\n";
  }
  epoch_times.clear();
}

void clear_buffer_stats() {
  for(unsigned int i = 0; i < flushes_size; ++i) 
    flushes[i] = 0;
  full = 0;
  messages_received = 0;
}

void clear_cumulative_buffer_stats() {
  for(unsigned int i = 0; i < flushes_size; ++i) 
    cumulative_flushes[i] = 0;
  cumulative_full = 0;
  cumulative_messages = 0;
}

void print_buffer_stats() {
  unsigned long long sum = 0;
  unsigned long long number_of_buffers = 0;

  flush_type tmp_flushes[flushes_size];
  for(int i = 0; i < flushes_size; ++i) 
    tmp_flushes[i] = flushes[i].load();
  flush_type tmp_full = full.load();
  flush_type tmp_messages = messages_received.load();

  MPI_Allreduce(&tmp_flushes, &all_flushes.front(), flushes_size, FLUSH_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&tmp_full, &all_full, 1, FLUSH_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&tmp_messages, &all_messages, 1, FLUSH_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);

  if(_RANK == 0) std::cout << "Flushes: ";
  for(unsigned int i = 0; i < all_flushes.size(); ++i) {
    if(all_flushes[i] != 0) {
      if(_RANK == 0) std::cout << i+1 << "->" << all_flushes[i] << "  ";
      sum += (i + 1) * all_flushes[i];
      number_of_buffers += all_flushes[i];
      cumulative_flushes[i] += all_flushes[i];
    }
  }
  cumulative_full += all_full;
  cumulative_messages += all_messages;
  if(_RANK == 0) {
    std::cout << std::endl;
    std::cout << "There were " << number_of_buffers << " incomplete flushes with the average size of a flush of " << (number_of_buffers != 0 ? sum/number_of_buffers : 0) << std::endl;
    std::cout << "Full buffers: " << all_full << std::endl;
    std::cout << "Messages: " << all_messages << std::endl;
  }
}

namespace amplusplus {
  namespace performance_counters {
    void hook_flushed_message_size(amplusplus::rank_type dest, size_t count, size_t elt_size)
    {
      if(count <= flushes_size)
	flushes[count-1]++;

#ifdef PRINT_STATS
      time_type t = get_time();
      fprintf(stderr, "Flush: %d bytes to %d at %f\n", (count * elt_size), dest, t);
#endif
    }

    void hook_full_buffer_send(amplusplus::rank_type dest, size_t count, size_t elt_size)
    {
      full++;
#ifdef PRINT_STATS
      time_type t = get_time();
      fprintf(stderr, "Full buffer: %d bytes to %d at %f\n", (count * elt_size), dest, t);
#endif
    }

    void hook_message_received(amplusplus::rank_type src, size_t count, size_t elt_size)
    {
      messages_received += count;
#ifdef PRINT_STATS
      time_type t = get_time();
      fprintf(stderr, "%d: Message received: %d bytes from %d at %f\n", _RANK, (count * elt_size), src, t);
#endif
    }

    void hook_begin_epoch(amplusplus::transport& trans)
    {
      epoch_times.push_back(get_time());
    }

    void hook_epoch_finished(amplusplus::transport& trans)
    {
      epoch_times.back() = get_time() - epoch_times.back();
    }

  }
}
#endif // AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS

// Hit rates
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
unsigned long long cumulative_hits, cumulative_tests;

void print_and_accumulate_cache_stats(std::pair<unsigned long long, unsigned long long> stats) {
  unsigned long long hits, tests;
  MPI_Allreduce(&stats.first, &hits, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stats.second, &tests, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  cumulative_hits += hits;
  cumulative_tests += tests;
  if(_RANK == 0) {
    std::cout << "Cache hit rates: " << hits << " hits, " << tests << " tests. Success ratio: " << (double)hits / (double)tests << "." << std::endl;
  }
}
#endif // AMPLUSPLUS_PRINT_HIT_RATES

#ifdef PBGL2_PRINT_WORK_STATS
typedef std::tuple<unsigned long long, unsigned long long, unsigned long long, unsigned long long> work_stats_t;

// Distributed control work stats
work_stats_t dc_stats, ds_stats;

void clear_cumulative_work_stats() {
  dc_stats = work_stats_t{ 0ul, 0ul, 0ul, 0ul };
  ds_stats = work_stats_t{ 0ul, 0ul, 0ul, 0ul };
}

void print_and_accumulate_work_stats(work_stats_t stats, work_stats_t &cummulative_stats) {
  work_stats_t temp_stats;
  MPI_Allreduce(&std::get<0>(stats), &std::get<0>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&std::get<1>(stats), &std::get<1>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&std::get<2>(stats), &std::get<2>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&std::get<3>(stats), &std::get<3>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(_RANK == 0) std::cout << "Useful work: " << std::get<0>(temp_stats) << ", ivalidated work: " << std::get<1>(temp_stats) << ", useless work: " << std::get<2>(temp_stats) << ", rejected work: " << std::get<3>(temp_stats) << std::endl;
  std::get<0>(cummulative_stats) += std::get<0>(temp_stats);
  std::get<1>(cummulative_stats) += std::get<1>(temp_stats);
  std::get<2>(cummulative_stats) += std::get<2>(temp_stats);
  std::get<3>(cummulative_stats) += std::get<3>(temp_stats);
}
#define PBGL2_DC_PRINT \
<< "\nUseful work: " << (std::get<0>(dc_stats) / num_sources) << " (per source), ivalidated work: " << (std::get<1>(dc_stats) / num_sources) << " (per source), useless work: " << (std::get<2>(dc_stats)/num_sources) << " (per source), rejected work: " << (std::get<3>(dc_stats)/num_sources) << " (per source). Rejected/useful ratio: " << ((double)std::get<3>(dc_stats) / (double)std::get<0>(dc_stats)) << ". Invalidated/useful ratio: " << ((double)std::get<1>(dc_stats) / (double)std::get<0>(dc_stats))
#define PBGL2_DS_PRINT \
<< "\nUseful work: " << (std::get<0>(ds_stats) / num_sources) << " (per source), ivalidated work: " << (std::get<1>(ds_stats) / num_sources) << " (per source), useless work: " << (std::get<2>(ds_stats)/num_sources) << " (per source), rejected work: " << (std::get<3>(ds_stats)/num_sources) << " (per source). Rejected/useful ratio: " << ((double)std::get<3>(ds_stats) / (double)std::get<0>(ds_stats)) << ". Invalidated/useful ratio: " << ((double)std::get<1>(ds_stats) / (double)std::get<0>(ds_stats))
#endif // PBGL2_PRINT_WORK_STATS

//
// Edge properties
//
typedef int32_t weight_type; // Needed so atomics are available on MacOS

struct WeightedEdge {
  WeightedEdge(weight_type weight = 0) : weight(weight) { }
  
  weight_type weight;
};

namespace amplusplus {
  template<>
  struct make_mpi_datatype<WeightedEdge> : make_mpi_datatype_base {
    make_mpi_datatype<weight_type> dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1() {
      int blocklengths[1] = {1};
      MPI_Aint displacements[1];
      WeightedEdge test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.weight, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Datatype types[1] = {dt1.get()};
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}


//
// Edge weight generator iterator
//
template<typename F, typename RandomGenerator>
class generator_iterator
{
public:
  typedef std::input_iterator_tag iterator_category;
  typedef typename F::result_type value_type;
  typedef const value_type&       reference;
  typedef const value_type*       pointer;
  typedef void                    difference_type;

  explicit 
  generator_iterator(RandomGenerator& gen, const F& f = F()) 
    : f(f), gen(&gen) 
  { 
    value = this->f(gen); 
  }

  reference operator*() const  { return value; }
  pointer   operator->() const { return &value; }

  generator_iterator& operator++()
  {
    value = f(*gen);
    return *this;
  }

  generator_iterator operator++(int)
  {
    generator_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const generator_iterator& other) const
  { return f == other.f; }

  bool operator!=(const generator_iterator& other) const
  { return !(*this == other); }

private:
  F f;
  RandomGenerator* gen;
  value_type value;
};

template<typename F, typename RandomGenerator>
inline generator_iterator<F, RandomGenerator> 
make_generator_iterator( RandomGenerator& gen, const F& f)
{ return generator_iterator<F, RandomGenerator>(gen, f); }

//
//
//
template <typename Graph, typename OwnerMap, typename LocalMap>
class cc_vertex_compare
{
  typedef typename property_traits<OwnerMap>::key_type Vertex;
  
public:
  cc_vertex_compare (const OwnerMap& o, const LocalMap& l)
    : owner(o), local(l) { }
  
  bool operator() (const Vertex& x, const Vertex& y) 
  { 
    if (x == graph_traits<Graph>::null_vertex()) return false;
    if (y == graph_traits<Graph>::null_vertex()) return true;

    if (get(local, x) < get(local, y))
      return true;
    else if (get(local, x) == get(local, y))
      return (get(owner, x) < get(owner, y));
    return false;
  }
  
private:
  OwnerMap   owner;
  LocalMap   local;
};


//
// BFS test helper
//
template <typename Graph, typename ColorMap, typename Queue>
time_type
run_bfs(amplusplus::transport& trans, amplusplus::transport& barrier_trans, 
	Graph& g, ColorMap& color, shared_ptr<Queue> Q,
	typename graph_traits<Graph>::vertex_descriptor current_source, int num_threads, 
	typename graph_traits<Graph>::vertices_size_type n, bool verify)
{
  typedef typename property_traits<ColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;

  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

  set_property_map_role(vertex_color, color);
  color.set_consistency_model(0);

  BGL_FORALL_VERTICES_T(v, g, Graph) { put(color, v, Color::white()); }

  trans.set_nthreads(num_threads);
  boost::graph::distributed::breadth_first_search<Graph, bfs_visitor<>, ColorMap, Queue> 
    bfs(g, make_bfs_visitor(null_visitor()), color, Q);
  bfs.set_source(current_source);
  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
  // num_threads threads now
  trans.set_nthreads(num_threads);
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif

  time_type start = get_time();
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(bfs), i + 1);
    threads[i].swap(thr);
  }

  // Run algorithm 
  bfs(0);

  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif // AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS

  // Back to one thread
  trans.set_nthreads(1);

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(color, v) == Color::black()) 
      ++visited; 
  }
  
  vertices_size_type total = boost::parallel::all_reduce_sum(barrier_trans, visited);

  if (verify) 
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << std::endl;

  // NGE: This is exceedingly arbitrary... remove later
  if (total <= 100) return -1.;

  return end - start;
}

struct flip_pair {
  template <typename T>
  T operator()(const T& x) const {
    return std::make_pair(x.second, x.first);
  }
};

template <typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_delta_stepping(amplusplus::transport& trans, amplusplus::transport& barrier_trans, Graph& g, const WeightMap& weight, 
		   DistanceMap& distance, typename graph_traits<Graph>::vertex_descriptor current_source, 
		   weight_type delta, int num_threads, typename graph_traits<Graph>::vertices_size_type n, bool verify, 
		   bool level_sync, MessageGenerator msg_gen) {
#ifdef CLONE
    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
    BGL_FORALL_VERTICES_T(v, g, Graph) 
      { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
    trans.set_nthreads(num_threads);
    boost::graph::distributed::delta_stepping_shortest_paths<Graph, DistanceMap, WeightMap,
	                                                     append_buffer<Vertex, 10u>, MessageGenerator>
      D(g, distance, weight, delta, msg_gen);
    trans.set_nthreads(1);

    D.set_source(current_source);
    
    if (level_sync)
      D.set_level_sync();

    { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
    time_type start = get_time();
    
    // Many threads now
    trans.set_nthreads(num_threads);
	  
    boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
    for (int i = 0; i < num_threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < num_threads - 1; ++i)
      threads[i].join();
	  
    time_type end = get_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    print_and_clear_epoch_times();
    print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
    const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
    print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
    work_stats_t work_stats = D.get_work_stats();
    print_and_accumulate_work_stats(work_stats, ds_stats);
#endif

    // Back to one thread
    trans.set_nthreads(1);

    unsigned long long num_levels = D.get_num_levels();

    if (verify) {

      //
      // TODO: DistanceMap uses trans which assumes num_threads are present! 
      //
	
      distance.set_consistency_model(boost::parallel::cm_forward);
      distance.set_max_ghost_cells(0);
	    
      {
	amplusplus::scoped_epoch epoch(g.transport());
	  
	BGL_FORALL_VERTICES_T(v, g, Graph) {
	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    get(distance, target(e, g));
	  }
	}
      }
	    
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  #ifdef PRINT_DEBUG
	  if (get(distance, target(e, g)) > boost::closed_plus<weight_type>()(get(distance, source(e, g)), get(weight, e)))
	    std::cout << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g)) << "->"
		      << get(get(vertex_local, g), target(e, g)) << "@" << get(get(vertex_owner, g), target(e, g)) << "  weight = "
		      << get(weight, e)
		      << "  distance(" << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g))
		      << ") = " << get(distance, source(e, g)) << "  distance(" << get(get(vertex_local, g), target(e, g)) << "@" 
		      << get(get(vertex_owner, g), target(e, g)) << ") = " << get(distance, target(e, g)) << std::endl;
	  #else
	  if(get(distance, target(e, g)) > 
	     boost::closed_plus<weight_type>()(get(distance, v), get(weight, e))) std::abort();
	  #endif
	}
      }
      
      if (trans.rank() == 0) std::cout << "Verified." << std::endl;
	    
      distance.clear(); // Clear memory used by ghost cells
    }	    

    vertices_size_type visited = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (get(distance, v) < std::numeric_limits<weight_type>::max())
	++visited; 
    }
	  
    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
      r(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r(visited);
	  
    //if (verify)
      if (trans.rank() == 0)
	std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		  << ",  " << num_levels << " levels required." << std::endl;

    if (total < 100) return -1.;
	  
    return end - start;
}



template <typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator, typename PriorityQueueGenerator = boost::graph::distributed::default_priority_queue_gen>
time_type
run_distributed_control(amplusplus::transport& trans, amplusplus::transport& barrier_trans, Graph& g, const WeightMap& weight, 
			DistanceMap& distance,  typename graph_traits<Graph>::vertex_descriptor current_source, int num_threads,
			typename graph_traits<Graph>::vertices_size_type n, bool verify, MessageGenerator msg_gen, MessageGenerator priority_msg_gen, int flushFreq) {
#ifdef CLONE
    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

	// Initialize with infinite (max) distance
    BGL_FORALL_VERTICES_T(v, g, Graph) 
	{ put(distance, v, std::numeric_limits<weight_type>::max()); }

    // find the incoming edge factor - i.e. an approximation of number of incoming edges per vertex
    unsigned int in_edges =  (num_edges(g) / num_vertices(g)) / 2;

    trans.set_nthreads(num_threads);
    boost::graph::distributed::distributed_control<Graph, DistanceMap, WeightMap,  PriorityQueueGenerator, MessageGenerator>
      D(g, distance, weight, g.transport(), msg_gen, priority_msg_gen,flushFreq);
    trans.set_nthreads(1);
    //trans.set_recvdepth(recvDepth);

    D.set_source(current_source);
    

    { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
	epoch_times.clear();
	clear_buffer_stats();
#endif
    
    time_type start = get_time();
    
    // Many threads now
    trans.set_nthreads(num_threads);
	  
    boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
    for (int i = 0; i < num_threads - 1; ++i) {
		boost::thread thr(boost::ref(D), i + 1);
		threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < num_threads - 1; ++i)
		threads[i].join();
	  
    time_type end = get_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    print_and_clear_epoch_times();
    print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
    const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
    print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
    work_stats_t work_stats = D.get_work_stats();
    print_and_accumulate_work_stats(work_stats, dc_stats);
#endif

    // Back to one thread
    trans.set_nthreads(1);



    if (verify) {

		//
		// TODO: DistanceMap uses trans which assumes num_threads are present! 
		//
	
		distance.set_consistency_model(boost::parallel::cm_forward);
		distance.set_max_ghost_cells(0);
	    
		{
			amplusplus::scoped_epoch epoch(g.transport());
	  
			BGL_FORALL_VERTICES_T(v, g, Graph) {
				BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
					get(distance, target(e, g));
				}
			}
		}
	    
		BGL_FORALL_VERTICES_T(v, g, Graph) {
			BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
#ifdef PRINT_DEBUG
				if (get(distance, target(e, g)) > boost::closed_plus<weight_type>()(get(distance, source(e, g)), get(weight, e)))
					std::cout << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g)) << "->"
							  << get(get(vertex_local, g), target(e, g)) << "@" << get(get(vertex_owner, g), target(e, g)) << "  weight = "
							  << get(weight, e)
							  << "  distance(" << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g))
							  << ") = " << get(distance, source(e, g)) << "  distance(" << get(get(vertex_local, g), target(e, g)) << "@" 
							  << get(get(vertex_owner, g), target(e, g)) << ") = " << get(distance, target(e, g)) << std::endl;
#else
				if(get(distance, target(e, g)) > 
				   boost::closed_plus<weight_type>()(get(distance, source(e, g)), get(weight, e))) std::abort();
#endif
			}
		}
	
		if (trans.rank() == 0) std::cout << "Verified" << std::endl;
    
		distance.clear(); // Clear memory used by ghost cells
    }	    

    vertices_size_type visited = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
		if (get(distance, v) < std::numeric_limits<weight_type>::max())
			++visited; 
    }
	  
    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
		r(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r(visited);
	  
    //if (verify)
    if (trans.rank() == 0) {
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start)  << "\n" << std::endl;
    }

    if (total < 100) return -1.;
	
    // if (trans.rank() == 0) {
      // D.print_in_edge_map();
    // }

    return end - start;
}



template <typename Graph, typename WeightMap, typename DistanceMap, typename PredecessorMap, typename MessageGenerator>
time_type
run_self_stabilizing(amplusplus::transport& trans, amplusplus::transport& barrier_trans, Graph& g, 
                     const WeightMap& weight, DistanceMap& distance, PredecessorMap& predecessor, 
                     typename graph_traits<Graph>::vertex_descriptor current_source, 
		     int num_threads, typename graph_traits<Graph>::vertices_size_type n, bool verify, 
		     bool stats, MessageGenerator msg_gen, unsigned int ordering_size, unsigned int levels, unsigned int delta) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

#ifdef PRINT_DEBUG
  std::cerr << "run_self_stabilizing: started\n" << std::flush;
#endif
  
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    put(distance, v, std::numeric_limits<weight_type>::max()); 
    put(predecessor, v, v); // it does not really matter how we initialize that map as long as we can tell which values have been touched by the algorithm 
  }

#ifdef PRINT_DEBUG
  std::cerr << "run_self_stabilizing: setting the number of threads to " << num_threads << "\n" << std::flush;
#endif
  trans.set_nthreads(num_threads);
#ifdef PRINT_DEBUG
  std::cerr << "run_self_stabilizing: The number of threads has been successfully set to " << num_threads << "\n" << std::flush;
#endif
  boost::graph::distributed::self_stabilizing_shortest_paths<Graph, DistanceMap, WeightMap, PredecessorMap, MessageGenerator>
    S(g, distance, weight, predecessor, ordering_size, levels, delta, msg_gen);
  trans.set_nthreads(1);
  
  S.set_source(current_source);
  
  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_buffer_stats();
#endif

  time_type start = get_time();
  
  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(S), i + 1);
    threads[i].swap(thr);
  }

  S(0);

  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();

  time_type end = get_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_buffer_stats();
#endif

#ifdef PRINT_DEBUG    
  std::cerr << "We have joined all the threads\n" << std::flush;
#endif

  // Back to one thread
  trans.set_nthreads(1);
  
#ifdef PRINT_DEBUG
  std::cerr << "run_self_stabilizing: verify is " << verify << "\n" << std::flush;
#endif

  if (verify) {
  
    //
    // TODO: DistanceMap uses trans which assumes num_threads are present! 
    //
  
    distance.set_consistency_model(boost::parallel::cm_forward);
    distance.set_max_ghost_cells(0);
  
    {
      amplusplus::scoped_epoch epoch(g.transport());
  
      BGL_FORALL_VERTICES_T(v, g, Graph) {
        BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
          get(distance, target(e, g));
        }
      }
    }
	    
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	//#ifdef PRINT_DEBUG
	  if (get(distance, target(e, g)) > boost::closed_plus<weight_type>()(get(distance, source(e, g)), get(weight, e)))
	    std::cout << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g)) << "->"
		      << get(get(vertex_local, g), target(e, g)) << "@" << get(get(vertex_owner, g), target(e, g)) << "  weight = "
		      << get(weight, e)
		      << "  distance(" << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g))
		      << ") = " << get(distance, source(e, g)) << "  distance(" << get(get(vertex_local, g), target(e, g)) << "@" 
		      << get(get(vertex_owner, g), target(e, g)) << ") = " << get(distance, target(e, g)) << std::endl;
#if 0
	  //#else
	  assert(get(distance, target(e, g)) <= 
		 boost::closed_plus<weight_type>()(get(distance, source(e, g)), get(weight, e)));
#endif
      }
    }
	    
    distance.clear(); // Clear memory used by ghost cells
  }

  vertices_size_type visited = 0;
  if(stats) {
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (get(distance, v) < std::numeric_limits<weight_type>::max())
	++visited; 
    }
  }  
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r_count(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r_count(visited);
	  
  if (stats)
    if (trans.rank() == 0) {
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) << "." << std::endl;
    }

  if (total < 100) return -1.;
	  
  return end - start;
}
	    
template <typename Graph, typename ParentMap, typename ComponentMap, typename LockMap, typename MessageGenerator>
time_type
run_sv_cc(amplusplus::transport& trans, Graph& g, ParentMap& parent, ComponentMap& component,
	  LockMap& locks, int num_threads, bool verify, bool level_sync, MessageGenerator msg_gen) {

#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  BGL_FORALL_VERTICES_T(v, g, Graph) {
      put(parent, v, v);
  }

  trans.set_nthreads(num_threads);
  boost::graph::distributed::connected_components<Graph, ParentMap, MessageGenerator>
    CC(g, parent, locks, msg_gen);
  trans.set_nthreads(1);

  if (level_sync) 
    CC.set_level_sync();
	  
  { amplusplus::scoped_epoch epoch(trans); }
  
  time_type start = get_time();
	  
  // Many threads now
  trans.set_nthreads(num_threads);
	  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(CC), i + 1);
    threads[i].swap(thr);
  }
	  
  CC(0);
	  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
	  
  time_type end = get_time();
  
  // Back to one thread
  trans.set_nthreads(1);

#ifdef NUMBER_COMPONENTS
  int num_components = CC.template number_components<ComponentMap>(component);
  if (trans.rank() == 0)
    std::cout << num_components << " components found\n";
#endif
  
  if (verify) {
      
    parent.set_max_ghost_cells(0);
	      
    {
      amplusplus::scoped_epoch epoch(trans);
	      
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(parent, source(e, g));
	  get(parent, target(e, g));
	}
      }
    }

    {	    
      amplusplus::scoped_epoch epoch(trans); // at the moment get() sends a message

      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_ADJ_T(v, u, g, Graph) {
#ifdef PRINT_DEBUG
	  if (get(parent, v) != get(parent, u)) 
	    std::cout << trans.rank() << ": parent(" << get(get(vertex_local, g), v) << "@"
		      << get(get(vertex_owner, g), v) << ") = " << get(get(vertex_local, g), get(parent, v)) 
		      << "@" << get(get(vertex_owner, g), get(parent, v)) << "  parent("
		      << get(get(vertex_local, g), u) << "@" << get(get(vertex_owner, g), u) 
		      << ") = " << get(get(vertex_local, g), get(parent, u)) 
		      << "@" << get(get(vertex_owner, g), get(parent, u)) << std::endl;
#else
	  assert(get(parent, v) == get(parent, u)); 
#endif
	}
      }
    }
    
    parent.clear(); // Clear memory used by ghost cells
  }
	
  return time_type(end - start);
}

template <typename Graph, typename ParentMap, typename ComponentMap, typename LockMap, typename MessageGenerator>
time_type
run_ps_sv_cc(amplusplus::transport& trans, Graph& g, ParentMap& parent, ComponentMap& component,
	     LockMap& locks, int num_threads, bool verify, bool level_sync, MessageGenerator msg_gen) {

#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(parent, v, v); // Initialize to self because we're going to use null_vertex() as a sentinel value in parallel search
  }

  // Instantiate algorithms used later here so we don't time their constructions
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
	  
  using boost::parallel::all_reduce;
  using boost::parallel::maximum;

  all_reduce<vertices_size_type, maximum<vertices_size_type> > 
    reduce_max(trans, maximum<vertices_size_type>());

  // TODO: If level_sync we should make the PS below a BFS

  boost::graph::distributed::parallel_search<Graph, ParentMap, MessageGenerator> 
    PS(g, parent, graph_traits<Graph>::null_vertex(), msg_gen);

  trans.set_nthreads(num_threads);
  boost::graph::distributed::connected_components<Graph, ParentMap, MessageGenerator>
    CC(g, parent, locks, msg_gen);
  trans.set_nthreads(1);

  if (level_sync) 
    CC.set_level_sync();

  
  // Less than on vertices w/ special casing for null_vertex() so we start parallel 
  // search from the same vertex below
  typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  typedef typename property_map<Graph, vertex_local_t>::const_type LocalMap;
  typedef cc_vertex_compare<Graph, OwnerMap, LocalMap> VertexLessThan;
  VertexLessThan vertex_lt(get(vertex_owner, g), get(vertex_local, g)); 

  //
  // Start timing loop
  //
  { amplusplus::scoped_epoch epoch(trans); }
	
  time_type start = get_time();

  //
  // Find max-degree vertex (gotta time this part too)
  //
  Vertex max_v = graph_traits<Graph>::null_vertex(); // Avoid uninitialized warning
  {
    vertices_size_type degree = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (out_degree(v, g) > degree) {
	degree = out_degree(v, g);
	max_v = v;
      }
    }
	  
    using boost::parallel::all_reduce;
    using boost::parallel::maximum;
	  
    vertices_size_type max = reduce_max(degree);
    
    // If the max-degree vertex is somewhere else then clear local max v
    // Note: This fails to handle two vertices with the same degree... in 
    //       that case the parallel search may start from more than one 
    //       vertex, which is perfectly acceptable provided they're in the
    //       same component... a reasonable but not fool-proof assumption
    if (max != degree) 
      max_v = graph_traits<Graph>::null_vertex();
  }

  // Find global minimum vertex w/ max degree as source	so all procs start 
  // at the same place
  all_reduce<Vertex, VertexLessThan> min_vertex(trans, vertex_lt);
  max_v = min_vertex(max_v);

  //
  // Run parallel search from max_v, this will mark all reachable vertices with 
  // null_vertex() in parent map
  //
  PS.run(max_v);

#ifdef PRINT_STATS
  all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    reduce_plus(trans, std::plus<vertices_size_type>());

  vertices_size_type giant_component_count = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (get(parent, v) == graph_traits<Graph>::null_vertex())
      ++giant_component_count;
  }
  vertices_size_type total = reduce_plus(giant_component_count);

  if (trans.rank() == 0)	
    std::cout << total << " vertices in giant component\n";
#endif

  //
  // Run SV CC on filtered graph
  //
  trans.set_nthreads(num_threads);
	
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(CC), i + 1);
    threads[i].swap(thr);
  }
	
  CC(0);
	
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
	
  time_type end = get_time();
	
  // Back to one thread
  trans.set_nthreads(1);
  
#ifdef NUMBER_COMPONENTS
  int num_components = CC.template number_components<ComponentMap>(component);
  if (trans.rank() == 0)
    std::cout << num_components << " components found\n";
#endif

  if (verify) {
      
    parent.set_max_ghost_cells(0);
    
    {
      amplusplus::scoped_epoch epoch(trans);
      
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(parent, source(e, g));
	  get(parent, target(e, g));
	}
      }
    }
	  
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_ADJ_T(v, u, g, Graph) {
#ifdef PRINT_DEBUG
	if (get(parent, v) != get(parent, u))
	  std::cout << trans.rank() << ": parent(" << get(get(vertex_local, g), v) << "@"
		    << get(get(vertex_owner, g), v) << ") = " << get(get(vertex_local, g), get(parent, v)) 
		    << "@" << get(get(vertex_owner, g), get(parent, v)) << "  parent("
		    << get(get(vertex_local, g), u) << "@" << get(get(vertex_owner, g), u) 
		    << ") = " << get(get(vertex_local, g), get(parent, u)) 
		    << "@" << get(get(vertex_owner, g), get(parent, u)) << std::endl;
#else
	assert(get(parent, v) == get(parent, u)); 
#endif
      }
    }
    
    parent.clear(); // Clear memory used by ghost cells
  }

  return time_type(end - start);
}

template <typename Graph, typename RankMap, typename MessageGenerator>
time_type
run_page_rank(amplusplus::transport& trans, Graph& g, RankMap& rank, RankMap& rank2, size_t iterations, 
	      float damping, typename graph_traits<Graph>::vertices_size_type n, int num_threads) 
{
  trans.set_nthreads(num_threads);
  boost::graph::distributed::page_rank<Graph, RankMap, boost::graph::n_iterations>
    PR(g, rank, rank2, boost::graph::n_iterations(iterations), damping, n);
  
  { amplusplus::scoped_epoch epoch(trans); }
  trans.set_nthreads(1);  

  time_type start = get_time();
  
  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(PR), i + 1);
    threads[i].swap(thr);
  }
	  
  PR(0);
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();
  
  // Back to one thread
  trans.set_nthreads(1);

  return end - start;
}

enum mode_type {mode_none, mode_self, mode_async_bfs, mode_ds_async_bfs, mode_level_synchronized_bfs, mode_delta_stepping, mode_dc, mode_connected_components, mode_page_rank};
enum routing_type {rt_none, rt_hypercube, rt_rook};
 
//
// main
//
int main(int argc, char* argv[]) {
//   {
//     int i = 0;
//     char hostname[256];
//     gethostname(hostname, sizeof(hostname));
//     printf("PID %d on %s ready for attach\n", getpid(), hostname);
//     fflush(stdout);
//     while (0 == i)
//       sleep(5);
//   }
 size_t recvDepth=1;
for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if(arg == "--receiveDepth") {
      recvDepth = boost::lexical_cast<size_t>( argv[i+1] );
      //std::cerr<<"receiveDepth from main"<<recvDepth<<std::endl;
    }
    
 }
 amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true, recvDepth);
  amplusplus::transport trans = env.create_transport();
  amplusplus::transport barrier_trans = trans.clone();

  _RANK = trans.rank(); // For performance counter output
	
  std::vector<int> thread_num_vals;
  size_t scale = 20, num_sources = 64, iterations = 20;
  double edgefactor = 16;
  unsigned long long n = static_cast<unsigned long long>(floor(pow(2, scale)));
  bool verify = false, level_sync = false, stats = true;
  uint64_t seed64 = 12345;
  weight_type C = 100;
  weight_type delta_min = 1, delta_max = 1, delta_step = 1;
  unsigned int levels_min = 1, levels_max = 1, levels_step = 1;
  mode_type mode = mode_none;
  routing_type routing = rt_none;
  double edge_list_reserve_factor = 1.15;
  bool no_reductions = false, per_thread_reductions = true;
  size_t coalescing_size = 4096, distribution_coalescing_size = 1 << 17;
  size_t reduction_cache_size = 20;
  size_t vertices_per_lock = 64;
  unsigned int ordering_size = 4098;
  unsigned int number_poll_task=1;
  std::string distributedControlDataStructure;
  int flushFreq=2000;
  size_t priority_coalescing_size = 40000;
 
  
#ifdef PRINT_ET
  time_type job_start = get_time();
#endif

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "--threads") {
      size_t d = 0, d2;
      std::string s = argv[i+1];

      while ((d2 = s.find(',', d)) != std::string::npos) {
	thread_num_vals.push_back(boost::lexical_cast<int>(s.substr(d, d2 - d)));
	d = d2 + 1;
      }
	
      thread_num_vals.push_back(boost::lexical_cast<int>(s.substr(d, s.length())));
    }

    if (arg == "--seed") {
	seed64 = boost::lexical_cast<uint64_t>( argv[i+1] );
    }

    if (arg == "--scale") {
      scale = boost::lexical_cast<size_t>( argv[i+1] );
      n = (unsigned long long)(1) << scale;
    }

    if (arg == "--degree")
      edgefactor = boost::lexical_cast<double>( argv[i+1] );

    if (arg == "--num-sources")
      num_sources = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--iterations")
      iterations = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--vertices-per-lock")
      vertices_per_lock = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--verify")
      verify = true;

    if (arg == "--level-sync")
      level_sync = true;

    if (arg == "--stats")
      stats = true;

    if (arg == "--coalescing-size") {
      coalescing_size = boost::lexical_cast<size_t>( argv[i+1] );
    }

    if (arg == "--ordering-size") {
      ordering_size = boost::lexical_cast<unsigned int>( argv[i+1] );
    }

    if (arg == "--reduction-cache-size")
        reduction_cache_size = boost::lexical_cast<size_t>( argv[i+1] );
    
    if (arg == "--distribution-coalescing-size")
	distribution_coalescing_size = boost::lexical_cast<size_t>( argv[i+1] );

    if (arg == "--max-weight")
      C = boost::lexical_cast<weight_type>( argv[i+1] );
    
    if(arg == "--poll_task"){
	    number_poll_task = boost::lexical_cast<unsigned int>( argv[i+1] );
	    //std::cerr<<"poll task number:" <<number_poll_task<<"\n";
    }

    if (arg == "--delta") {
      std::string s = argv[i+1];
      size_t d = s.find(':');
      if(d == std::string::npos) { std::cout << "Invalid value of delta\n"; abort(); }
      size_t d2 = s.find(':', d+1);
      if(d2 == std::string::npos) { std::cout << "Invalid value of delta\n"; abort(); }

      delta_min = boost::lexical_cast<weight_type>(s.substr(0, d));
      delta_max = boost::lexical_cast<weight_type>(s.substr(d+1, d2 - d - 1));
      delta_step = boost::lexical_cast<weight_type>(s.substr(d2 + 1, s.length()));

      if (delta_min > delta_max) { std::cout << "Invalid value of delta\n"; abort(); }
      if (delta_step <= 0) { std::cout << "Invalid value of delta\n"; abort(); }
    }
    
     if (arg == "--flush") {
	     flushFreq =  boost::lexical_cast<int> ( argv[i+1] );
     }
     
     if (arg == "--levels") {
      std::string s = argv[i+1];
      size_t d = s.find(':');
      if(d == std::string::npos) { std::cout << "Invalid value of levels\n"; abort(); }
      size_t d2 = s.find(':', d+1);
      if(d2 == std::string::npos) { std::cout << "Invalid value of levels\n"; abort(); }

      levels_min = boost::lexical_cast<weight_type>(s.substr(0, d));
      levels_max = boost::lexical_cast<weight_type>(s.substr(d+1, d2 - d - 1));
      levels_step = boost::lexical_cast<weight_type>(s.substr(d2 + 1, s.length()));

      if (levels_min > levels_max) { std::cout << "Invalid value of levels\n"; abort(); }
      if (levels_step <= 0) { std::cout << "Invalid value of levels\n"; abort(); }
    }

    if (arg == "--self") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_self;
    }

    if (arg == "--bfs") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_level_synchronized_bfs;
    }

    if (arg == "--async-bfs") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_async_bfs;
    }

    if (arg == "--ds-async-bfs") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_ds_async_bfs;
    }

    if (arg == "--delta-stepping") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_delta_stepping;
    }

	if (arg == "--dc" || arg == "--distributed_control") { // distributed_control
		if (mode != mode_none) 
			goto multiple_mode_settings; 
		mode = mode_dc;
    }

	if(arg == "--ds"){
		distributedControlDataStructure = argv[i+1];
	}
	       
	
		

    if (arg == "--connected-components") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_connected_components;
    }

    if (arg == "--page-rank") {
      if (mode != mode_none) 
	goto multiple_mode_settings; 
      mode = mode_page_rank;
    }

    if (arg == "--hypercube")
      routing = rt_hypercube;

    if (arg == "--rook")
      routing = rt_rook;

    if (arg == "--reserve-factor")
      edge_list_reserve_factor = boost::lexical_cast<double>( argv[i+1] );

    if (arg == "--with-no-reductions")
      no_reductions = true;

    if (arg == "--without-no-reductions")
      no_reductions = false;

    if (arg == "--with-per-thread-reductions")
      per_thread_reductions = true;

    if (arg == "--without-per-thread-reductions")
      per_thread_reductions = false;
  }


#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  flushes = new atomic_flush_type[coalescing_size]{};
  flushes_size = coalescing_size;
  all_flushes.resize(coalescing_size);
  cumulative_flushes.resize(coalescing_size);
#endif

  if (mode == mode_none) {
      std::clog << (boost::format("%d: need to set a mode of either --self, --dc or --distributed_control, --bfs, --async-bfs, --ds-async-bfs, --delta-stepping, --page-rank, or --connected-components\n") % argv[0]).str() << std::flush;
    return 1;

    if (0) {
      multiple_mode_settings:
        std::clog << (boost::format("%d: need to set exactly one mode\n") % argv[0]).str() << std::flush;
      return 1;
    }
  }

  if (thread_num_vals.empty()) thread_num_vals.push_back(1);
#if 0
  bool threads_required = false;
  for (std::vector<int>::iterator iter = thread_num_vals.begin() ; iter != thread_num_vals.end() ; ++iter)
      if (*iter != 1)
	  threads_required = true;

  // AM++ initialization                                                                                                                                                                                                                       
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, threads_required, recvDepth, number_poll_task);
  //std::cerr<<"receiveDepth after environment initialization"<<recvDepth<<std::endl;
  amplusplus::transport trans = env.create_transport();
  amplusplus::transport barrier_trans = trans.clone();
#endif
  typedef amplusplus::transport::rank_type rank_type;

  // Seed general-purpose RNG
  rand48 gen, synch_gen;
  gen.seed(seed64);
  synch_gen.seed(seed64);

  //
  // Directed graph for BFS variants
  //
  if (mode == mode_level_synchronized_bfs || mode == mode_async_bfs || mode == mode_ds_async_bfs || mode == mode_page_rank) {

    typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property, 
      distributedS<unsigned long long> >
      Digraph;

    typedef graph_traits<Digraph>::vertex_descriptor Vertex;
    typedef graph_traits<Digraph>::edges_size_type edges_size_type;
    typedef graph_traits<Digraph>::vertices_size_type vertices_size_type;

    typedef property_map<Digraph, vertex_index_t>::type VertexIndexMap;
    typedef property_map<Digraph, vertex_owner_t>::const_type OwnerMap;

    edges_size_type m = static_cast<edges_size_type>(floor(n * edgefactor));

    time_type gen_start = get_time();

    // Output of permutation
    std::vector<std::pair<edges_size_type, edges_size_type> > edges;
    edges.reserve(static_cast<edges_size_type>(floor(edge_list_reserve_factor * 2 * m / trans.size())));

    parallel::variant_distribution<vertices_size_type> distrib 
      = parallel::block<vertices_size_type>(trans, n);

    {
#ifdef CLONE
      amplusplus::transport trans = trans.clone(); // Clone transport for distribution
#endif
      // Build a graph to test with
      typedef graph500_iterator<Digraph> Graph500Iter;

      boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());
      
      edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
      edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

      // Permute and redistribute copy constructs the input iterator
      uint64_t a = rand_64(gen);
      uint64_t b = rand_64(gen);

      //
      // Distribute graph using specified routing policy
      //
      if (routing == rt_none) {

        typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::no_routing(trans.rank(), trans.size()));

	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	
      if (routing == rt_hypercube) {

        typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	  
      if (routing == rt_rook) {

        typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
    }

#ifdef BGP_REPORT_MEMORY
    long long memSize, usedMem; 

#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
    getMemSize(&memSize);
    getUsedMem(&usedMem);
    std::cout << trans.rank() << ": Edge distribution complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

    Digraph g(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), n, trans, distrib);

    // Clear edge array above
    edges.clear();

    // Generate sources
    boost::uniform_int<uint64_t> rand_vertex(0, n-1);
      
    { amplusplus::scoped_epoch epoch(trans); }

    time_type gen_end = get_time();

    if (trans.rank() == 0)
      std::cout << "Graph generation took " << print_time(gen_end - gen_start) << "s\n\n";

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
    getMemSize(&memSize);
    getUsedMem(&usedMem);
    std::cout << trans.rank() << ": Graph construction complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

    //
    // Scope for level-synchronized BFS
    //
    if (mode == mode_level_synchronized_bfs) {

      typedef boost::parallel::global_index_map<VertexIndexMap,
	property_map<Digraph, vertex_global_t>::type>
	GlobalIndexMap;
      GlobalIndexMap global_index(trans, num_vertices(g), get(vertex_index, g), get(vertex_global, g));
    
#ifdef PRINT_MAX_DEGREE_VERTEX
      {
	typedef graph_traits<Digraph>::vertices_size_type vertices_size_type;
      
	Vertex max_v = vertex(0, g);
	vertices_size_type degree = 0;
	BGL_FORALL_VERTICES(v, g, Digraph) { 
	  if (out_degree(v, g) > degree) {
	    degree = out_degree(v, g);
	    max_v = v;
	  }
	}
	
	using boost::parallel::all_reduce;
	using boost::parallel::maximum;
	
	all_reduce<edges_size_type, maximum<edges_size_type> > r(trans, maximum<edges_size_type>());
	edges_size_type max = r(degree);
	
	if (max == degree)
	  std::cout << trans.rank() << ": Max degree vertex (index = " << get(global_index, max_v) 
		    << ", degree = " << degree << ")\n";
      }
#endif

      for (unsigned int i = 0; i < thread_num_vals.size(); ++i) {
        const unsigned int num_threads = thread_num_vals[i];

	trans.set_nthreads(num_threads); // DEBUG
	typedef two_bit_color_map<VertexIndexMap> ColorMap;
	ColorMap color(num_vertices(g), get(vertex_index, g));
	trans.set_nthreads(1); // DEBUG

	set_property_map_role(vertex_color, color);
	color.set_consistency_model(0);

	time_type total_bfs_time = 0;
	std::pair<unsigned long long, unsigned long long> cache_stats;

	if (per_thread_reductions) {
	  //
	  // Call Level-synchronized BFS with no global color caching
	  //
#ifdef BFS_VECTOR
	  if (num_threads == 1) {

	    if (routing == rt_none) {
#ifdef CLONE
	      amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	      
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats()
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
		  std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		using boost::graph::distributed::distributed_queue;
		typedef distributed_queue<OwnerMap, std::vector<Vertex>,
	                                  boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                  MessageGenerator>
		  distributed_queue_type;
	  
		trans.set_nthreads(num_threads);
		shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
		  					  		        boost::make_shared<std::vector<Vertex> >(),
									        boost::make_shared<std::vector<Vertex> >(),
									        boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									        msg_gen));
		trans.set_nthreads(1);

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
		time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

		if (bfs_time == -1.) { // Not enough vertices visited
		  --source_i; continue;
		} 

#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
		std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
		cache_stats.first += t.first;
		cache_stats.second += t.second;
#endif

		Q.reset();
		
		total_bfs_time += bfs_time;
	      }
	    }
	  
	    if (routing == rt_hypercube) {
#ifdef CLONE
	      amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
		
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats()
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
		  std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		using boost::graph::distributed::distributed_queue;
		typedef distributed_queue<OwnerMap, std::vector<Vertex>,
	                                  boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                MessageGenerator>
		  distributed_queue_type;
	  
		trans.set_nthreads(num_threads);
		shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
		  	  						        boost::make_shared<std::vector<Vertex> >(),
									        boost::make_shared<std::vector<Vertex> >(),
									        boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									        msg_gen));
		trans.set_nthreads(1);

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
		
		time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      
	      

		if (bfs_time == -1.) { // Not enough vertices visited
		  --source_i; continue;
		} 

#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
		std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
		cache_stats.first += t.first;
		cache_stats.second += t.second;
#endif

		Q.reset();
		
		total_bfs_time += bfs_time;
	      }
	    }

	    if (routing == rt_rook) {
#ifdef CLONE		
	      amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	      
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats()
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
		  std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		using boost::graph::distributed::distributed_queue;
		typedef distributed_queue<OwnerMap, std::vector<Vertex>,
	                                  boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                  MessageGenerator>
		  distributed_queue_type;
	  
		trans.set_nthreads(num_threads);
		shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
		  							        boost::make_shared<std::vector<Vertex> >(),
									        boost::make_shared<std::vector<Vertex> >(),
									        boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									        msg_gen));
		trans.set_nthreads(1);

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
		
		time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      
	      
		if (bfs_time == -1.) { // Not enough vertices visited
		  --source_i; continue;
		} 

#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
		std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
		cache_stats.first += t.first;
		cache_stats.second += t.second;
#endif

		Q.reset();
	      
		total_bfs_time += bfs_time;
	      }
	    }

	    if (trans.rank() == 0)
	      std::cout << "Total BFS (no global caching, duplicate elimination, std::vector) time for " << num_sources << " sources = " 
			<< print_time(total_bfs_time) << " (" 
			<< print_time(total_bfs_time / num_sources) << " per source)  " << num_threads << " threads\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
			<< "\nTests = " << cache_stats.first << ", hits = " << cache_stats.second
			<< ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n";
#else
	    ;
#endif
	  } // end if (per_thread_reductions)

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 
          getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif

#endif // BFS_VECTOR

  	  //
	  // Call Level-synchronized BFS with no global color caching
	  //
	  if (per_thread_reductions) {
	    total_bfs_time = 0;
	    cache_stats = std::make_pair(0,0);

	    if (routing == rt_none) {
#ifdef CLONE
	      amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
		  std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		using boost::graph::distributed::distributed_queue;
		typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                  boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                  MessageGenerator>
	          distributed_queue_type;

		trans.set_nthreads(num_threads);
		shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  						        make_shared<append_buffer<Vertex, 20> >(),
									        make_shared<append_buffer<Vertex, 20> >(),
									        boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									        msg_gen));
		trans.set_nthreads(1);

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
		time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

		if (bfs_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

#if 0
		std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
		cache_stats.first += t.first;
		cache_stats.second += t.second;
#endif

		Q.reset();
		
		total_bfs_time += bfs_time;
	      }
	    }

	    if (routing == rt_hypercube) {
#ifdef CLONE
	      amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	      
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
		  std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		using boost::graph::distributed::distributed_queue;
		typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                  boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                  MessageGenerator>
	          distributed_queue_type;

		trans.set_nthreads(num_threads);
		shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			    						        make_shared<append_buffer<Vertex, 20> >(),
									        make_shared<append_buffer<Vertex, 20> >(),
									        boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									        msg_gen));
		trans.set_nthreads(1);

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
		time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

		if (bfs_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 
		
#if 0
		std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
		cache_stats.first += t.first;
		cache_stats.second += t.second;
#endif
		
		Q.reset();
		
		total_bfs_time += bfs_time;
	      }
	    }

	    if (routing == rt_rook) {
#ifdef CLONE
	      amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	      
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
		  std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		using boost::graph::distributed::distributed_queue;
		typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                  boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                  MessageGenerator>
	          distributed_queue_type;

		trans.set_nthreads(num_threads);
		shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			    						        make_shared<append_buffer<Vertex, 20> >(),
									        make_shared<append_buffer<Vertex, 20> >(),
									        boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									        msg_gen));
		trans.set_nthreads(1);
		
		Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
		time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

		if (bfs_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

#if 0
		std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
		cache_stats.first += t.first;
		cache_stats.second += t.second;
#endif
		
		Q.reset();
		
		total_bfs_time += bfs_time;
	      }
	    } 
	    
	    if (trans.rank() == 0)
	      std::cout << "Total BFS (no global caching, duplicate elimination, append_buffer) time for " << num_sources << " sources = " 
			<< print_time(total_bfs_time) << " (" 
			<< print_time(total_bfs_time / num_sources) << " per source)  " << num_threads << " threads\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
			<< "\nTests = " << cache_stats.first << ", hits = " << cache_stats.second
			<< ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n";
#else
	      ;
#endif
	  }

#ifdef BGP_REPORT_MEMORYY
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 
          getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif
	} // end if (per_thread_reductions)

	//
	// Call Level-synchronized BFS with no global color caching or reductions
	//
	if (no_reductions) {
	  total_bfs_time = 0;
	  cache_stats = std::make_pair(0,0);

	  if (routing == rt_none) {
#ifdef CLONE
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	    MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {
   	    
#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	
	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                MessageGenerator>
	        distributed_queue_type;

	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  					  	      make_shared<append_buffer<Vertex, 20> >(),
									      make_shared<append_buffer<Vertex, 20> >(),
									      boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									      msg_gen));

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

#if 0
	      std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
	      cache_stats.first += t.first;
	      cache_stats.second += t.second;
#endif

	      Q.reset();

	      total_bfs_time += bfs_time;
	    }
	  }

	  if (routing == rt_hypercube) {
#ifdef CLONE		
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {
   	    
#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	
	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                MessageGenerator>
	        distributed_queue_type;

	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  					  	      make_shared<append_buffer<Vertex, 20> >(),
									      make_shared<append_buffer<Vertex, 20> >(),
									      boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									     msg_gen));

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

#if 0
	      std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
	      cache_stats.first += t.first;
	      cache_stats.second += t.second;
#endif

	      Q.reset();

	      total_bfs_time += bfs_time;
	    }
	  }

	  if (routing == rt_rook) {
#ifdef CLONE		
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {
   	    
#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	
	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>,
	                                MessageGenerator>
	        distributed_queue_type;

	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  					  	      make_shared<append_buffer<Vertex, 20> >(),
									      make_shared<append_buffer<Vertex, 20> >(),
									      boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), trans.rank()),
									     msg_gen));

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

#if 0
	      std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
	      cache_stats.first += t.first;
	      cache_stats.second += t.second;
#endif

	      Q.reset();

	      total_bfs_time += bfs_time;
	    }
	  }

	  if (trans.rank() == 0)
	    std::cout << "Total BFS (no global caching, no reductions) time for " << num_sources << " sources = " 
		      << print_time(total_bfs_time) << " (" 
		      << print_time(total_bfs_time / num_sources) << " per source)  " << num_threads << " threads\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "\nTests = " << cache_stats.first << ", hits = " << cache_stats.second
		      << ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n";
#else
	      ;
#endif
	}

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
        if (trans.rank() == 0) {
#endif 
        getMemSize(&memSize);
	getUsedMem(&usedMem);
	std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 
#endif

#ifndef DISABLE_BFS_BITMAP
	if (per_thread_reductions) {
	  //
	  // Call Level-synchronized BFS with global bitmap
	  //
	  total_bfs_time = 0;
	  cache_stats = std::make_pair(0,0);

	  if (routing == rt_none) {
#ifdef CLONE	
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex>,
       	                                boost::detail::has_not_been_seen<GlobalIndexMap>,
 	                                MessageGenerator>
	        distributed_queue_type;
	  
	      trans.set_nthreads(num_threads);
	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  						      make_shared<append_buffer<Vertex> >(),
									      make_shared<append_buffer<Vertex> >(),
									      boost::detail::has_not_been_seen<GlobalIndexMap>(n, global_index),// get(vertex_index, g)),
									      msg_gen));
	      trans.set_nthreads(1);

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 
	  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
	      std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
	      cache_stats.first += t.first;
	      cache_stats.second += t.second;
#endif
	  
	      Q.reset();
	      
	      total_bfs_time += bfs_time;
	    }
	  }

	  if (routing == rt_hypercube) {
#ifdef CLONE	
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex>,
       	                                boost::detail::has_not_been_seen<GlobalIndexMap>,
 	                                MessageGenerator>
	        distributed_queue_type;
	  
	      trans.set_nthreads(num_threads);
	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  						      make_shared<append_buffer<Vertex> >(),
									      make_shared<append_buffer<Vertex> >(),
									      boost::detail::has_not_been_seen<GlobalIndexMap>(n, global_index),// get(vertex_index, g)),
									      msg_gen));
	      trans.set_nthreads(1);

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 
	  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
	      std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
	      cache_stats.first += t.first;
	      cache_stats.second += t.second;
#endif
	  
	      Q.reset();
	      
	      total_bfs_time += bfs_time;
	    }
	  }

	  if (routing == rt_rook) {
#ifdef CLONE	
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex>,
       	                                boost::detail::has_not_been_seen<GlobalIndexMap>,
 	                                MessageGenerator>
	        distributed_queue_type;
	  
	      trans.set_nthreads(num_threads);
	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
			  						      make_shared<append_buffer<Vertex> >(),
									      make_shared<append_buffer<Vertex> >(),
									      boost::detail::has_not_been_seen<GlobalIndexMap>(n, global_index),// get(vertex_index, g)),
									      msg_gen));
	      trans.set_nthreads(1);

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 
	  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
	      std::pair<unsigned long long, unsigned long long> t = Q->get_cache_stats();
	      cache_stats.first += t.first;
	      cache_stats.second += t.second;
#endif
	  
	      Q.reset();
	      
	      total_bfs_time += bfs_time;
	    }
	  }
	
	  if (trans.rank() == 0)
	    std::cout << "Total BFS (global bitmap, duplicate elimination) time for " << num_sources << " sources = " 
		      << print_time(total_bfs_time) << " (" 
		      << print_time(total_bfs_time / num_sources) << " per source)  " << num_threads << " threads\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "\nTests = " << cache_stats.first << ", hits = " << cache_stats.second
		      << ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n";
#else
	  ;
#endif
	} // end if (per_thread_reductions)

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
        if (trans.rank() == 0) {
#endif 
        getMemSize(&memSize);
	getUsedMem(&usedMem);
	std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 
#endif

	if (no_reductions) {
	  //
	  // Call Level-synchronized BFS with global bitmap
	  //
	  total_bfs_time = 0;
	  cache_stats = std::make_pair(0,0);

	  if (routing == rt_none) {
#ifdef CLONE
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                boost::detail::has_not_been_seen<GlobalIndexMap>, MessageGenerator>
	        distributed_queue_type;
	  
	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
									      make_shared<append_buffer<Vertex, 20> >(),
									      make_shared<append_buffer<Vertex, 20> >(),
									      boost::detail::has_not_been_seen<GlobalIndexMap>(n, global_index),// get(vertex_index, g)),
									      msg_gen));
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 
	  
	      Q.reset();

	      total_bfs_time += bfs_time;
	    }
	  }

	  if (routing == rt_hypercube) {
#ifdef CLONE
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                boost::detail::has_not_been_seen<GlobalIndexMap>, MessageGenerator>
	        distributed_queue_type;
	  
	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
									      make_shared<append_buffer<Vertex, 20> >(),
									      make_shared<append_buffer<Vertex, 20> >(),
									      boost::detail::has_not_been_seen<GlobalIndexMap>(n, global_index),// get(vertex_index, g)),
									      msg_gen));
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 
	  
	      Q.reset();

	      total_bfs_time += bfs_time;
	    }
	  }

	  if (routing == rt_rook) {
#ifdef CLONE
	    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      using boost::graph::distributed::distributed_queue;
	      typedef distributed_queue<OwnerMap, append_buffer<Vertex, 20>,
	                                boost::detail::has_not_been_seen<GlobalIndexMap>, MessageGenerator>
	        distributed_queue_type;
	  
	      shared_ptr<distributed_queue_type> Q(new distributed_queue_type(trans, get(vertex_owner, g), 
									      make_shared<append_buffer<Vertex, 20> >(),
									      make_shared<append_buffer<Vertex, 20> >(),
									      boost::detail::has_not_been_seen<GlobalIndexMap>(n, global_index),// get(vertex_index, g)),
									      msg_gen));
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type bfs_time = run_bfs<Digraph, ColorMap, distributed_queue_type>
#ifdef BARRIER_TRANS
		  (trans, barrier_trans, g, color, Q, current_source, num_threads, n, verify);
#else
		  (trans, trans, g, color, Q, current_source, num_threads, n, verify);
#endif	      

	      if (bfs_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 
	  
	      Q.reset();

	      total_bfs_time += bfs_time;
	    }
	  }
	
	  if (trans.rank() == 0)
	    std::cout << "Total BFS (global bitmap, no reductions) time for " << num_sources << " sources = " 
		      << print_time(total_bfs_time) << " (" 
		      << print_time(total_bfs_time / num_sources) << " per source)  " << num_threads << " threads\n";
	}

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
        if (trans.rank() == 0) {
#endif 
        getMemSize(&memSize);
	getUsedMem(&usedMem);
	std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 
#endif

#endif // DISABLE_BFS_BITMAP
      }
    }

    //
    // Asynchronous BFS
    //
    if (mode == mode_async_bfs) {
      // Distance map
      std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
      typedef iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap> 
	DistanceMap;
      DistanceMap distance(distanceS.begin(), get(vertex_index, g));


      {
	//
	// Asynchronous BFS no reductions
	//

	time_type total_async_bfs_time = 0;
	std::pair<unsigned long long, unsigned long long> msg_stats(0,0);
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {
	  
	  BGL_FORALL_VERTICES(v, g, Digraph)
	    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
	  typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	  MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	  
	  boost::graph::distributed::async_breadth_first_search<Digraph, DistanceMap, MessageGenerator>
	    async_bfs(g, distance, msg_gen);
	  
	  { amplusplus::scoped_epoch epoch(trans); }
	  
	  time_type start = get_time();
	  
	  Vertex current_source = vertex(rand_vertex(synch_gen), g);

	  async_bfs.run(current_source);
	  
	  time_type end = get_time();

	  total_async_bfs_time += end - start;

#ifdef AMPLUSPLUS_PRINT_HIT_RATES	  
	  std::pair<unsigned long long, unsigned long long> t = async_bfs.get_msg_stats();
	  msg_stats.first += t.first;
	  msg_stats.second += t.second;
#endif

	  if (verify) {
	    vertices_size_type visited = 0;
	    BGL_FORALL_VERTICES(v, g, Digraph) { 
	      if (get(distance, v) < std::numeric_limits<weight_type>::max())
		++visited; 
	    }
	    
	    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
	      r(trans, std::plus<vertices_size_type>());
	    vertices_size_type total = r(visited);
	    
	    if (trans.rank() == 0)
	      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) << std::endl;
	  }
	}
	
	if (trans.rank() == 0)
	  std::cout << "Total Async BFS (no reductions) time for " << num_sources << " sources = " 
		    << print_time(total_async_bfs_time) << " (" 
		    << print_time(total_async_bfs_time / num_sources) << " per source)\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		    << "  MESSAGE Tests = " << msg_stats.first << ", hits = " << msg_stats.second
		    << ", hit rate = " << (double(msg_stats.second) / double(msg_stats.first) * 100.) << "%\n";
#else
	;
#endif    
      
#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	if (trans.rank() == 0) {
#endif 
	getMemSize(&memSize);
	getUsedMem(&usedMem);
	std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 
#endif
	total_async_bfs_time = 0;
	std::pair<unsigned long long, unsigned long long> cache_stats(0,0);

	if (per_thread_reductions) {

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
   	  //
	  // Asynchronous BFS w/ reductions
	  //
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {
	  
	    BGL_FORALL_VERTICES(v, g, Digraph)
	      { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	  
	    boost::graph::distributed::async_breadth_first_search<Digraph, DistanceMap, MessageGenerator>
	      async_bfs(g, distance, msg_gen);
	  
	    { amplusplus::scoped_epoch epoch(trans); }
	    
	    time_type start = get_time();
	    
	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    
	    async_bfs.run(current_source);
	  
	    time_type end = get_time();
	    
	    total_async_bfs_time += end - start;
	  
#if 0
	    std::pair<unsigned long long, unsigned long long> t = async_bfs.get_cache_stats();
	    cache_stats.first += t.first;
	    cache_stats.second += t.second;
	    
	    t = async_bfs.get_msg_stats();
	    msg_stats.first += t.first;
	    msg_stats.second += t.second;
#endif
	  
	    if (verify) {
	      vertices_size_type visited = 0;
	      BGL_FORALL_VERTICES(v, g, Digraph) { 
	        if (get(distance, v) < std::numeric_limits<weight_type>::max())
		  ++visited; 
	      }
	    
	      boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
	        r(trans, std::plus<vertices_size_type>());
	      vertices_size_type total = r(visited);
	    
	      if (trans.rank() == 0)
	        std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) << std::endl;
	    }
	  }
	
	  if (trans.rank() == 0)
	    std::cout << "Total Async BFS (w/ reductions) time for " << num_sources << " sources = " 
		      << print_time(total_async_bfs_time) << " (" 
		      << print_time(total_async_bfs_time / num_sources) << " per source)\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "  CACHE Tests = " << cache_stats.first << ", hits = " << cache_stats.second
		      << ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n"
		      << "  MESSAGE Tests = " << msg_stats.first << ", hits = " << msg_stats.second
		      << ", hit rate = " << (double(msg_stats.second) / double(msg_stats.first) * 100.) << "%\n";
#else
	  ;
#endif

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 
          getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif
	} // end if (per_thread_reductions)
      }
    }

    //
    // Delta-Stepping Asynchronous BFS
    //
    if (mode == mode_ds_async_bfs) {

      // Distance map
      std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
      typedef iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap> 
	DistanceMap;
      DistanceMap distance(distanceS.begin(), get(vertex_index, g));

      // TODO: Unit weight map
      typedef static_property_map<weight_type> WeightMap;
      WeightMap weight(static_cast<weight_type>(1));

      for (unsigned int i = 0, num_threads = thread_num_vals[i]; i < thread_num_vals.size(); ++i, num_threads = thread_num_vals[i]) {
	for (weight_type delta = delta_min; delta <= delta_max; delta += delta_step) {

	  if (trans.rank() == 0) 
	    std::cout << trans.rank() << ": Delta = " << delta << std::endl;

	  time_type total_delta_stepping_time = 0;
	  std::pair<unsigned long long, unsigned long long> cache_stats(0,0);

	  if (no_reductions) {
	    //
	    // Delta-stepping BFS no reductions
	    //

	    total_delta_stepping_time = 0;
	    cache_stats = std::make_pair(0, 0);

	    if (routing == rt_none) {

	      typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	      MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
	          std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
		
		time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	          (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

		if (ds_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

		total_delta_stepping_time += ds_time;
	      }
	    }

	    if (routing == rt_hypercube) {

	      typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
	          std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
		time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	          (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

		if (ds_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

		total_delta_stepping_time += ds_time;
	      }
	    }

	    if (routing == rt_rook) {

	      typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
	          std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		Vertex current_source = vertex(rand_vertex(synch_gen), g);

		time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	          (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

		if (ds_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

		total_delta_stepping_time += ds_time;
	      }
	    }

	    if (trans.rank() == 0)
	      std::cout << "Total Delta-Stepping BFS (no reductions) time for " << num_sources << " sources = " 
		        << print_time(total_delta_stepping_time) << " (" 
		        << print_time(total_delta_stepping_time / num_sources) << " per source)  delta = " << delta 
		        << ", " << num_threads << " threads, " << "coalescing_size = " << coalescing_size << std::endl;

	  }

	  if (per_thread_reductions) {
	  //
	  // Delta-stepping BFS with reductions
	  //

	    total_delta_stepping_time = 0;
	    cache_stats = std::make_pair(0, 0);
	    
	    if (routing == rt_none) {

	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
	          std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
		Vertex current_source = vertex(rand_vertex(synch_gen), g);
		
		time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	          (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

		if (ds_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

		total_delta_stepping_time += ds_time;
	      }
	    }

	    if (routing == rt_hypercube) {

	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
	          std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

		Vertex current_source = vertex(rand_vertex(synch_gen), g);
		
		time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	          (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

		if (ds_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

		total_delta_stepping_time += ds_time;
	      }
	    }
	    
	    if (routing == rt_rook) {

	      typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	      MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	      for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	        if (trans.rank() == 0)
	          std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
		Vertex current_source = vertex(rand_vertex(synch_gen), g);
		
		time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	          (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

		if (ds_time == -1.) { // Not enough vertices visited
	          --source_i; continue;
		} 

		total_delta_stepping_time += ds_time;
	      }
	    }
	
	    if (trans.rank() == 0)
              std::cout << "Total Delta-Stepping BFS (w/ per-thread reductions) time for " << num_sources << " sources = " 
		        << print_time(total_delta_stepping_time) << " (" 
		        << print_time(total_delta_stepping_time / num_sources) << " per source)  delta = " << delta 
		        << ", " << num_threads << " threads, "<< "coalescing_size = " << coalescing_size << std::endl
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		        << "  Tests = " << cache_stats.first << ", hits = " << cache_stats.second
		        << ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n";
#else
	    ;
#endif
	  } // end if (per_thred_reductions)
	} // end delta loop
      } // end nthreads loop
    } // end mode_ds_async_bfs 

    if (mode == mode_page_rank) {

      typedef property_map<Digraph::base_type, vertex_index_t>::type LocalVertexIndexMap;
      typedef vector_property_map<float, LocalVertexIndexMap> LocalRankMap;
     
      typedef property_map<Digraph, vertex_global_t>::const_type VertexGlobalMap;

      // TODO: Only valid for 1 thread now
//       for (unsigned int i = 0, num_threads = thread_num_vals[i]; i < thread_num_vals.size(); ++i, num_threads = thread_num_vals[i]) {
      {
	unsigned int num_threads = 1;

        if (no_reductions) {
          //
          // Delta-stepping no reductions
          //

// 	  amplusplus::transport trans2 = trans.clone();

	  time_type pr_time = 0;
	
	  if (routing == rt_none) {

	    typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	    MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	
	    typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap, MessageGenerator> RankMap;
	    RankMap rank(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen),
	      rank2(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen);
	    
	    pr_time = run_page_rank<Digraph, RankMap, MessageGenerator>
	      (trans, g, rank, rank2, iterations, 0.85, n, num_threads);
	  }

	  if (routing == rt_hypercube) {

	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
	    typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap, MessageGenerator> RankMap;
	    RankMap rank(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen),
	      rank2(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen);

	    pr_time = run_page_rank<Digraph, RankMap, MessageGenerator>
	      (trans, g, rank, rank2, iterations, 0.85, n, num_threads);
	  }
	
	  if (routing == rt_rook) {

	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	
	    typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap, MessageGenerator> RankMap;
	    RankMap rank(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen),
	      rank2(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen);

	    pr_time = run_page_rank<Digraph, RankMap, MessageGenerator>
	      (trans, g, rank, rank2, iterations, 0.85, n, num_threads);
	  }

	  if (trans.rank() == 0)
	    std::cout << "Total PageRank (no reductions) time for " << iterations << " iterations = " 
		      << print_time(pr_time) << "   " << num_threads << " threads\n";
	}

	if (per_thread_reductions) {
	  //
          // Delta-stepping with reductions
          //

          time_type pr_time = 0;
	  std::pair<unsigned long long, unsigned long long> cache_stats(0,0);

	  if (routing == rt_none) {
	    
	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	    
	    typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap, MessageGenerator> RankMap;
	    RankMap rank(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen),
	      rank2(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen);

	    pr_time = run_page_rank<Digraph, RankMap, MessageGenerator>
	      (trans, g, rank, rank2, iterations, 0.85, n, num_threads);
	  }

	  if (routing == rt_hypercube) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	    typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap, MessageGenerator> RankMap;
	    RankMap rank(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen),
	      rank2(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen);

	    pr_time = run_page_rank<Digraph, RankMap, MessageGenerator>
	      (trans, g, rank, rank2, iterations, 0.85, n, num_threads);
	  }

	  if (routing == rt_rook) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	    
	    typedef boost::parallel::distributed_property_map<VertexGlobalMap, LocalRankMap, MessageGenerator> RankMap;
	    RankMap rank(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen),
	      rank2(trans, get(vertex_global, g), LocalRankMap(num_vertices(g), get(vertex_index, g.base())), msg_gen);

	    pr_time = run_page_rank<Digraph, RankMap, MessageGenerator>
	      (trans, g, rank, rank2, iterations, 0.85, n, num_threads);
	  }
	
	  if (trans.rank() == 0)
            std::cout << "Total PageRank (w/ per-thread reductions) time for " << iterations << " iterations = " 
		      << print_time(pr_time) << "   " << num_threads << " threads\n"
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "  Tests = " << cache_stats.first << ", hits = " << cache_stats.second
		      << ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%\n";
#else
	    ;
#endif
	} // end if (per_thread_reductions)
      } // end thread_num loop
    } // end mode_page_rank
  } // scope for graph used in BFSs


  //============================================================================//
  //
  // Distributed control
  //
  if (mode == mode_dc) {

    typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge, no_property, distributedS<unsigned long long> > Digraph;

    typedef graph_traits<Digraph>::vertex_descriptor Vertex;
    typedef graph_traits<Digraph>::edges_size_type edges_size_type;
    typedef graph_traits<Digraph>::vertices_size_type vertices_size_type;

    typedef property_map<Digraph, vertex_index_t>::type VertexIndexMap;
    typedef property_map<Digraph, vertex_owner_t>::const_type OwnerMap;

    edges_size_type m = static_cast<edges_size_type>(floor(n * edgefactor));

    time_type gen_start = get_time();

    // Output of permutation
    std::vector<std::pair<edges_size_type, edges_size_type> > edges;
    edges.reserve(static_cast<edges_size_type>(floor(edge_list_reserve_factor * 2 * m / trans.size())));

    parallel::variant_distribution<vertices_size_type> distrib = parallel::block<vertices_size_type>(trans, n);

    {
      // Build a graph to test with
      typedef graph500_iterator<Digraph> Graph500Iter;

      boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());

#ifdef CLONE
      amplusplus::transport trans = trans.clone(); // Clone transport for distribution
#endif
      
      edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
      edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

      // Permute and redistribute copy constructs the input iterator
      uint64_t a = rand_64(gen);
      uint64_t b = rand_64(gen);

      //
      // Distribute graph using specified routing policy
      //
      if (routing == rt_none) {

	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::no_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	
      if (routing == rt_hypercube) {

	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	  
      if (routing == rt_rook) {

	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
    }

#ifdef BGP_REPORT_MEMORY
    long long memSize, usedMem; 
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
      getMemSize(&memSize);
      getUsedMem(&usedMem);
      std::cout << trans.rank() << ": Edge distribution complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

    Digraph g(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), 
	      make_generator_iterator(gen, uniform_int<weight_type>(1, C)), n, trans, distrib);

    // Clear edge array above
    edges.clear();

    // Generate sources
    boost::uniform_int<uint64_t> rand_vertex(0, n-1);
      
    { amplusplus::scoped_epoch epoch(trans); }

    time_type gen_end = get_time();

    if (trans.rank() == 0)
      std::cout << "Graph generation took " << print_time(gen_end - gen_start) << "s\n";

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
      getMemSize(&memSize);
      getUsedMem(&usedMem);
      std::cout << trans.rank() << ": Graph construction complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

#ifdef PRINT_ET
    if (trans.rank() == 0)
      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

    // Property maps
    typedef property_map<Digraph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = get(&WeightedEdge::weight, g);
    
    // Distance map
    std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
    typedef iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap>  DistanceMap;
    DistanceMap distance(distanceS.begin(), get(vertex_index, g));

    // Processed edge vector - Used for priority
    /* std::vector<unsigned int> initial_processed_edges(num_vertices(g), 0); // initially processed edges are 0 for all vertices
    typedef iterator_property_map<std::vector<unsigned int>::iterator, VertexIndexMap> ProcessedInEdgeCountMap;
    ProcessedInEdgeCountMap processed_in_edge_counts(initial_processed_edges.begin(), get(vertex_index, g));
    */
    for (unsigned int i = 0, num_threads = thread_num_vals[i]; i < thread_num_vals.size(); ++i, num_threads = thread_num_vals[i]) {
	
      time_type total_dc_time = 0;

#ifdef DS_SHARED_REDUCTIONS
      //
      // Distributed control with reductions
      //

      // MZ: I am not sure if we are supposed to be hitting this code (see the #ifdef above).  As for now, make sure we are not.
      std::abort();

      if (routing == rt_none) {

	typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	  if (trans.rank() == 0)
	    std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	  Vertex current_source = vertex(rand_vertex(synch_gen), g);
	  time_type dc_time;
	  if(distributedControlDataStructure == "vector") {			 
	    dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
	      (trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	  }
	  else{
	    dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator>
	      (trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq );
		
	  }


	  if (dc_time == -1.) { // Not enough vertices visited
	    --source_i; continue;
	  } 

	  total_dc_time += dc_time;
	}
      }

      if (routing == rt_hypercube) {

	typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	  if (trans.rank() == 0)
	    std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	  Vertex current_source = vertex(rand_vertex(synch_gen), g);
	  time_type dc_time;
	  if(distributedControlDataStructure == "vector") {
	    dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
	      (trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	  }
	  else {
	    dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator>
	      (trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	  }



	  if (dc_time == -1.) { // Not enough vertices visited
	    --source_i; continue;
	  } 

	  total_dc_time += dc_time;
				  
	}
      }

      if (routing == rt_rook) {


	typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	  if (trans.rank() == 0)
	    std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	  Vertex current_source = vertex(rand_vertex(synch_gen), g);
	  time_type dc_time;
	  if(distributedControlDataStructure == "vector"){
	    dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
	      (trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	  }
	  else{
	    dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator>
	      (trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	  }

	  if (dc_time == -1.) { // Not enough vertices visited
	    --source_i; continue;
	  } 

	  total_dc_time += dc_time;
	}
      }

      if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	unsigned long long buffers{}, messages{};
	for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	  buffers += cumulative_flushes[i];
	  messages += cumulative_flushes[i] * i;
	}
#endif

	std::cout << "Total Distributed Control (DC) (w/ shared reductions) time for " << num_sources << " sources = " 
		  << print_time(total_dc_time) << " (" 
		  << print_time(total_dc_time / num_sources) << " per source), " << num_threads << " threads, coalescing_size =" << coalescing_size
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		  << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		  << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		  << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		  << (messages/buffers) << " messages per buffer on average."
#endif
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		  << "\nTests = " << (cumulative_tests / num_sources) << " (per source), hits = " << (cumulative_hits / num_sources)
		  << " (per source), hit rate = " << (double(cumulative_hits) / double(cumulative_tests) * 100.) << "%"
#endif
#ifdef PBGL2_PRINT_WORK_STATS
	  PBGL2_DC_PRINT // The print macro for work stats
#endif
		  << std::endl;
      }

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
      if (trans.rank() == 0) {
#endif 
	getMemSize(&memSize);
	getUsedMem(&usedMem);
	std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
      }
#endif 
#endif

#endif // DS_SHARED_REDUCTIONS

      if (no_reductions) {
	//
	// Distributed Control (DC) no reductions
	//

	total_dc_time = 0;

	if (routing == rt_none) {

	  typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	  MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	  MessageGenerator priority_msg_gen((amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1)));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    time_type dc_time; 
	    if(distributedControlDataStructure == "vector") {
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
		(trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    else{
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    if (dc_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_dc_time += dc_time;
					  
	  }
	}

	if (routing == rt_hypercube) {

	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    time_type dc_time ;
	    if(distributedControlDataStructure == "vector"){					  
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    else{
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }

	    if (dc_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_dc_time += dc_time;
					  
	  }
	}

	if (routing == rt_rook) {

	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    time_type dc_time ;
	    if(distributedControlDataStructure == "vector"){
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
		(trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    else{
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }


	    if (dc_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_dc_time += dc_time;

	  }
	}

	if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	  unsigned long long buffers{}, messages{};
	  for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	    buffers += cumulative_flushes[i];
	    messages += cumulative_flushes[i] * i;
	  }
#endif
	  std::cout << "Total Distributed Control (DC) (no reductions) time for " << num_sources << " sources = " 
		    << print_time(total_dc_time) << " (" 
		    << print_time(total_dc_time / num_sources) << " per source), " << num_threads << " threads, coalescing_size = " << coalescing_size 
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		    << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		    << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		    << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		    << (messages/buffers) << " messages per buffer on average."
#endif		      
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		    << "\nTests = " << (cumulative_tests / num_sources) << " (per source), hits = " << (cumulative_hits / num_sources)
		  << " (per source), hit rate = " << (double(cumulative_hits) / double(cumulative_tests) * 100.) << "%"
#endif
#ifdef PBGL2_PRINT_WORK_STATS
	  PBGL2_DC_PRINT // The print macro for work stats
#endif
		    << std::endl;
	}
#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	if (trans.rank() == 0) {
#endif 
	  getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 

#endif
      }

      if (per_thread_reductions) {
	//
	// Distributed Control with reductions
	//

	total_dc_time = 0;

	if (routing == rt_none) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    time_type dc_time ;
	    if(distributedControlDataStructure == "vector"){
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    else{
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator>
		(trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    if (dc_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_dc_time += dc_time;
	  }
	}

	if (routing == rt_hypercube) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));

	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    time_type dc_time ;
	    if(distributedControlDataStructure == "vector") {
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    else {
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }

	    if (dc_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_dc_time += dc_time;

	  }
	}

	if (routing == rt_rook) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));

	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";

#endif
	  
	    Vertex current_source = vertex(rand_vertex(synch_gen), g);
	    time_type dc_time ;
	    if(distributedControlDataStructure == "vector") {
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap,  MessageGenerator, boost::graph::distributed::vector_of_vector_gen>
		(trans, barrier_trans, g, weight, distance,  current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }
	    else {
	      dc_time = run_distributed_control<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq);
	    }

	    if (dc_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_dc_time += dc_time;

	  }
	}
	
	if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	  unsigned long long buffers{}, messages{};
	  for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	    buffers += cumulative_flushes[i];
	    messages += cumulative_flushes[i] * i;
	  }
#endif
	  float msg_per_buf = -1;
	  if (buffers != 0) {
	    msg_per_buf = messages/buffers;
	  }
				  
	  std::cout << "Total Distributed Control (DC) (w/ per-thread reductions) time for " << num_sources << " sources = " 
		    << print_time(total_dc_time) << " (" 
		    << print_time(total_dc_time / num_sources) << " per source), " << num_threads << " threads, coalescing_size = "
		    << coalescing_size << std::endl
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		    << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		    << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		    << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		    << msg_per_buf << " messages per buffer on average (-1 means buffer size is 0)."
#endif
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		    << "\nTests = " << (cumulative_tests / num_sources) << " (per source), hits = " << (cumulative_hits / num_sources)
		    << " (per source), hit rate = " << (double(cumulative_hits) / double(cumulative_tests) * 100.) << "%"
#endif
#ifdef PBGL2_PRINT_WORK_STATS
	    PBGL2_DC_PRINT // The print macro for work stats
#endif
		    << std::endl;

	}

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	if (trans.rank() == 0) {
#endif 
	  getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 
#endif
      } // end if (per_thread_reductions)

    } // end for thread_num loop
  }
  //
  // End of Distributed Control (DC)
  // 
  
  //
  // Delta-stepping shortest paths
  //
  if (mode == mode_delta_stepping) {

    typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge, no_property, 
					distributedS<unsigned long long> >
      Digraph;

    typedef graph_traits<Digraph>::vertex_descriptor Vertex;
    typedef graph_traits<Digraph>::edges_size_type edges_size_type;
    typedef graph_traits<Digraph>::vertices_size_type vertices_size_type;

    typedef property_map<Digraph, vertex_index_t>::type VertexIndexMap;
    typedef property_map<Digraph, vertex_owner_t>::const_type OwnerMap;

    edges_size_type m = static_cast<edges_size_type>(floor(n * edgefactor));

    time_type gen_start = get_time();

    // Output of permutation
    std::vector<std::pair<edges_size_type, edges_size_type> > edges;
    edges.reserve(static_cast<edges_size_type>(floor(edge_list_reserve_factor * 2 * m / trans.size())));

    parallel::variant_distribution<vertices_size_type> distrib = parallel::block<vertices_size_type>(trans, n);

    {
      // Build a graph to test with
      typedef graph500_iterator<Digraph> Graph500Iter;

      boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());

#ifdef CLONE
      amplusplus::transport trans = trans.clone(); // Clone transport for distribution
#endif
      
      edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
      edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

      // Permute and redistribute copy constructs the input iterator
      uint64_t a = rand_64(gen);
      uint64_t b = rand_64(gen);

      //
      // Distribute graph using specified routing policy
      //
      if (routing == rt_none) {

	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::no_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	
      if (routing == rt_hypercube) {

	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	  
      if (routing == rt_rook) {

	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
    }

#ifdef BGP_REPORT_MEMORY
    long long memSize, usedMem; 
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
      getMemSize(&memSize);
      getUsedMem(&usedMem);
      std::cout << trans.rank() << ": Edge distribution complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

    Digraph g(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), 
	      make_generator_iterator(gen, uniform_int<weight_type>(1, C)), n, trans, distrib);

    // Clear edge array above
    edges.clear();

    // Generate sources
    boost::uniform_int<uint64_t> rand_vertex(0, n-1);
      
    { amplusplus::scoped_epoch epoch(trans); }

    time_type gen_end = get_time();

    if (trans.rank() == 0)
      std::cout << "Graph generation took " << print_time(gen_end - gen_start) << "s\n";

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
      getMemSize(&memSize);
      getUsedMem(&usedMem);
      std::cout << trans.rank() << ": Graph construction complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

#ifdef PRINT_ET
    if (trans.rank() == 0)
      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

    // Property maps
    typedef property_map<Digraph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = get(&WeightedEdge::weight, g);
    
    // Distance map
    std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
    typedef iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap>  DistanceMap;
    DistanceMap distance(distanceS.begin(), get(vertex_index, g));

    for (unsigned int i = 0, num_threads = thread_num_vals[i]; i < thread_num_vals.size(); ++i, num_threads = thread_num_vals[i]) {
      for (weight_type delta = delta_min; delta <= delta_max; delta += delta_step) {
	
	time_type total_delta_stepping_time = 0;
	std::pair<unsigned long long, unsigned long long> cache_stats(0,0);

#ifdef DS_SHARED_REDUCTIONS
	//
	// Delta-stepping with reductions
	//

	if (routing == rt_none) {

	  typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	    Vertex current_source = vertex(rand_vertex(synch_gen), g);

	    time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	      (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	    if (ds_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_delta_stepping_time += ds_time;
	  }
	}

	if (routing == rt_hypercube) {

	  typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	    Vertex current_source = vertex(rand_vertex(synch_gen), g);

	    time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	      (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	    if (ds_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_delta_stepping_time += ds_time;
	  }
	}

	if (routing == rt_rook) {

	  typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	  for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	    if (trans.rank() == 0)
	      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	    Vertex current_source = vertex(rand_vertex(synch_gen), g);

	    time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
	      (trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	    if (ds_time == -1.) { // Not enough vertices visited
	      --source_i; continue;
	    } 

	    total_delta_stepping_time += ds_time;
	  }
	}

	if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS

	  unsigned long long buffers{}, messages{};
	  for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	    buffers += cumulative_flushes[i];
	    messages += cumulative_flushes[i] * i;
	  }
#endif


	  std::cout << "Total Delta-Stepping (w/ shared reductions) time for " << num_sources << " sources = " 
		    << print_time(total_delta_stepping_time) << " (" 
		    << print_time(total_delta_stepping_time / num_sources) << " per source)   delta = " << delta 
		    << ", " << num_threads << " threads, coalescing_size =" << coalescing_size

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		    << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		    << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		    << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		    << (messages/buffers) << " messages per buffer on average."
#endif
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		    << "\nTests = " << (cumulative_tests / num_sources) << " (per source), hits = " << (cumulative_hits / num_sources)
		  << " (per source), hit rate = " << (double(cumulative_hits) / double(cumulative_tests) * 100.) << "%"
#endif
#ifdef PBGL2_PRINT_WORK_STATS
	    PBGL2_DS_PRINT // The print macro for work stats
#endif
		    << std::endl;
}

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	if (trans.rank() == 0) {
#endif 
	  getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	}
#endif 
#endif

#endif // DS_SHARED_REDUCTIONS

	if (no_reductions) {
	  //
	  // Delta-stepping no reductions
	  //

	  total_delta_stepping_time = 0;
	  cache_stats = std::make_pair(0, 0);

	  if (routing == rt_none) {

	    typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	    MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	      if (ds_time == -1.) { // Not enough vertices visited
		--source_i; continue;
	      } 

	      total_delta_stepping_time += ds_time;
	    }
	  }

	  if (routing == rt_hypercube) {

	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	      if (ds_time == -1.) { // Not enough vertices visited
		--source_i; continue;
	      } 

	      total_delta_stepping_time += ds_time;
	    }
	  }

	  if (routing == rt_rook) {

	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	      if (ds_time == -1.) { // Not enough vertices visited
		--source_i; continue;
	      } 

	      total_delta_stepping_time += ds_time;
	    }
	  }

	  if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	    unsigned long long buffers{}, messages{};
	    for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	      buffers += cumulative_flushes[i];
	      messages += cumulative_flushes[i] * i;
	    }
#endif
	    std::cout << "Total Delta-Stepping (no reductions) time for " << num_sources << " sources = " 
		      << print_time(total_delta_stepping_time) << " (" 
		      << print_time(total_delta_stepping_time / num_sources) << " per source)  delta = " << delta 
		      << ", " << num_threads << " threads, coalescing_size = " << coalescing_size 
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		      << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		      << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		      << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		      << (messages/buffers) << " messages per buffer on average."
#endif		      
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "\nTests = " << (cumulative_tests / num_sources) << " (per source), hits = " << (cumulative_hits / num_sources)
		  << " (per source), hit rate = " << (double(cumulative_hits) / double(cumulative_tests) * 100.) << "%"
#endif
#ifdef PBGL2_PRINT_WORK_STATS
	    PBGL2_DS_PRINT // The print macro for work stats
#endif
		      << std::endl;
	  }

#ifdef BGP_REPORT_MEMORY
	  long long memSize, usedMem; 
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 

	    getMemSize(&memSize);
	    getUsedMem(&usedMem);
	    std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";

#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif
	}

	if (per_thread_reductions) {
	  //
	  // Delta-stepping with reductions
	  //

	  total_delta_stepping_time = 0;
	  cache_stats = std::make_pair(0, 0);

	  if (routing == rt_none) {


	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	      if (ds_time == -1.) { // Not enough vertices visited
		--source_i; continue;
	      } 

	      total_delta_stepping_time += ds_time;
	    }
	  }

	  if (routing == rt_hypercube) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	      if (ds_time == -1.) { // Not enough vertices visited
		--source_i; continue;
	      } 

	      total_delta_stepping_time += ds_time;
	    }
	  }

	  if (routing == rt_rook) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
		std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);
	      
	      time_type ds_time = run_delta_stepping<Digraph, WeightMap, DistanceMap, MessageGenerator>
		(trans, barrier_trans, g, weight, distance, current_source, delta, num_threads, n, verify, level_sync, msg_gen);

	      if (ds_time == -1.) { // Not enough vertices visited
		--source_i; continue;
	      } 

	      total_delta_stepping_time += ds_time;
	    }
	  }
	
	  if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	    unsigned long long buffers{}, messages{};
	    for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	      buffers += cumulative_flushes[i];
	      messages += cumulative_flushes[i] * i;
	    }
#endif
	    std::cout << "Total Delta-Stepping (w/ per-thread reductions) time for " << num_sources << " sources = " 
		      << print_time(total_delta_stepping_time) << " (" 
		      << print_time(total_delta_stepping_time / num_sources) << " per source)  delta = " << delta 
		      << ", " << num_threads << " threads, coalescing_size = " << coalescing_size
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		      << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		      << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		      << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		      << (messages/buffers) << " messages per buffer on average."
#endif
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "\nTests = " << (cumulative_tests / num_sources) << " (per source), hits = " << (cumulative_hits / num_sources)
		  << " (per source), hit rate = " << (double(cumulative_hits) / double(cumulative_tests) * 100.) << "%"
#endif
#ifdef PBGL2_PRINT_WORK_STATS
	    PBGL2_DS_PRINT // The print macro for work stats
#endif
		      << std::endl;
	  }

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 
	    getMemSize(&memSize);
	    getUsedMem(&usedMem);
	    std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif
	} // end if (per_thread_reductions)

      } // end for delta=delta_min:delta_max loop

    } // end for thread_num loop

  } // end all delta-stepping algorithms

  if (mode == mode_connected_components) {

    //
    // TODO: Create undirected ER Graph
    //
    typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property, 
                                        distributedS<unsigned long long> >
      Graph;

    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;

    typedef graph_traits<Graph>::edges_size_type edges_size_type;
    typedef graph_traits<Graph>::vertices_size_type vertices_size_type;

    time_type gen_start = get_time();

    // Output of permutation
    std::vector<std::pair<edges_size_type, edges_size_type> > edges, local_edges;
    local_edges.reserve(static_cast<size_t>(floor(edge_list_reserve_factor * 2 * edgefactor * n / trans.size())));
    edges.reserve(static_cast<size_t>(floor(edge_list_reserve_factor * 2 * edgefactor * n / trans.size())));

    parallel::variant_distribution<vertices_size_type> distrib 
      = parallel::block<vertices_size_type>(trans, n);

    // Build a graph to test with
    typedef erdos_renyi_iterator<rand48, Graph> ERIter;

    double p = (double)edgefactor / n;

    ERIter edge_iter(gen, n, p);
    std::advance(edge_iter, edgefactor * n * trans.rank() / trans.size()); // edges-per-proc * rank 

    for (unsigned int i = 0 ; i < edgefactor * n / trans.size() ; ++i) {
      std::pair<edges_size_type, edges_size_type> edge = *edge_iter;
      local_edges.push_back(edge);
      ++edge_iter;
    }

    //
    // Distribute graph using specified routing
    //
    {
#ifdef CLONE
      amplusplus::transport trans = trans.clone(); // Clone transport for distribution
#endif

      if (routing == rt_none) {

        typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::no_routing(trans.rank(), trans.size()));
	    
	distribute(local_edges.begin(), local_edges.end(), flip_pair(), std::back_inserter(edges), distrib, trans, msg_gen);
      }
	
      if (routing == rt_hypercube) {

        typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	distribute(local_edges.begin(), local_edges.end(), flip_pair(), std::back_inserter(edges), distrib, trans, msg_gen);
      }
	  
      if (routing == rt_rook) {

        typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	    
	distribute(local_edges.begin(), local_edges.end(), flip_pair(), std::back_inserter(edges), distrib, trans, msg_gen);
      }
    }

    // Clear local edges
    local_edges.clear();

#ifdef BGP_REPORT_MEMORY
    long long memSize, usedMem; 
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
    getMemSize(&memSize);
    getUsedMem(&usedMem);
    std::cout << trans.rank() << ": Edge distribution complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif
      
    Graph g(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), n, trans, distrib);

    // Clear edge array above
    edges.clear();

    // Generate sources
    boost::uniform_int<uint64_t> rand_vertex(0, n-1);
      
    { amplusplus::scoped_epoch epoch(trans); }

    time_type gen_end = get_time();

    if (trans.rank() == 0)
      std::cout << "Graph generation took " << print_time(gen_end - gen_start) << "s\n";

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
    getMemSize(&memSize);
    getUsedMem(&usedMem);
    std::cout << trans.rank() << ": Graph construction complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif
#endif

    // Parent, component, and lock maps for CC
    std::vector<Vertex> parentS(num_vertices(g), graph_traits<Graph>::null_vertex());
    typedef iterator_property_map<std::vector<Vertex>::iterator, VertexIndexMap> ParentMap;
    ParentMap parent(parentS.begin(), get(vertex_index, g));

#ifdef NUMBER_COMPONENTS
    std::vector<int> componentS(num_vertices(g), std::numeric_limits<int>::max());
    typedef iterator_property_map<std::vector<int>::iterator, VertexIndexMap> ComponentMap;
    ComponentMap component(componentS.begin(), get(vertex_index, g));
#endif

    typedef boost::parallel::lock_map<VertexIndexMap> LockMap;
    LockMap locks(get(vertex_index, g), num_vertices(g) / vertices_per_lock);
    
    time_type cc_time = 0;
    
    //
    // Connected components
    //
    for (unsigned int i = 0, num_threads = thread_num_vals[i]; i < thread_num_vals.size(); ++i, num_threads = thread_num_vals[i]) {

      //
      // Search in giant component, SV on filtered graph (all vertices not in giant component)
      //
      if (no_reductions) {
	    
        if (routing == rt_none) {
		  
	  typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	  MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	  
	  cc_time = run_ps_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	
	if (routing == rt_hypercube) {

  	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_ps_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	  
	if (routing == rt_rook) {
		
  	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_ps_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}

	if (trans.rank() == 0)
	  std::cout << "Total parallel search + SV CC (no reductions) time = " << cc_time << "  " << num_threads << " threads\n"; 
      }

      if (per_thread_reductions) {
	// Reductions in hooking and parallel search

        if (routing == rt_none) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	    
	  cc_time = run_ps_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	
	if (routing == rt_hypercube) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_ps_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	  
	if (routing == rt_rook) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_ps_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	    
	if (trans.rank() == 0)
	  std::cout << "Total parallel search + SV CC (per-thread reductions) time = " << cc_time << "  " << num_threads << " threads\n"; 
      }

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
      if (trans.rank() == 0) {
#endif 
	getMemSize(&memSize);
	getUsedMem(&usedMem);
	std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
      }
#endif
#endif

      //
      // Shiloach-Vishkin CC
      //

      if (no_reductions) {
	    
        if (routing == rt_none) {
		  
	  typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	  MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	  
	  cc_time = run_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}

	if (routing == rt_hypercube) {

  	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}

	if (routing == rt_rook) {
		
	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}

	if (trans.rank() == 0)
	  std::cout << "Total Shiloach-Vishkin CC (no reductions) time = " << cc_time << "  " << num_threads << " threads\n"; 
      }

      if (per_thread_reductions) {
	// Reductions in hooking

        if (routing == rt_none) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	  
	if (routing == rt_hypercube) {
	      
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	  cc_time = run_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	  
	if (routing == rt_rook) {

	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	  
	  cc_time = run_sv_cc(trans, g, parent, component, locks, num_threads, verify, level_sync, msg_gen);
	}
	    
	if (trans.rank() == 0)
	  std::cout << "Total Shiloach-Vishkin CC (per-thread reductions in hooking) time = " << cc_time << "  " << num_threads << " threads\n"; 
      }

    } // end for thread_num loop

  } // End CC mode loop

  if (mode == mode_self) {

#ifdef PRINT_DEBUG
    if(trans.rank() == 0)
      std::cout << "Mode is mode_self" << std::endl;
#endif

    typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge, no_property, 
                                        distributedS<unsigned long long> >
      Digraph;

    typedef graph_traits<Digraph>::vertex_descriptor Vertex;
    typedef graph_traits<Digraph>::edges_size_type edges_size_type;
    typedef graph_traits<Digraph>::vertices_size_type vertices_size_type;

    typedef property_map<Digraph, vertex_index_t>::type VertexIndexMap;
    typedef property_map<Digraph, vertex_owner_t>::const_type OwnerMap;

    edges_size_type m = static_cast<edges_size_type>(floor(n * edgefactor));

    time_type gen_start = get_time();

    // Output of permutation
    std::vector<std::pair<edges_size_type, edges_size_type> > edges;
    edges.reserve(static_cast<edges_size_type>(floor(edge_list_reserve_factor * 2 * m / trans.size())));

    parallel::variant_distribution<vertices_size_type> distrib = parallel::block<vertices_size_type>(trans, n);

#ifdef PRINT_DEBUG
    std::cerr << "europar: self stabilizing: about to generate the graph\n" << std::flush;
#endif

    {
      // Build a graph to test with
      typedef graph500_iterator<Digraph> Graph500Iter;

      boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());

#ifdef CLONE
      amplusplus::transport trans = trans.clone(); // Clone transport for distribution
#endif
      
      edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
      edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

      // Permute and redistribute copy constructs the input iterator
      uint64_t a = rand_64(gen);
      uint64_t b = rand_64(gen);

      //
      // Distribute graph using specified routing policy
      //
      if (routing == rt_none) {

#ifdef PRINT_DEBUG
	std::cerr << "europar: self stabilizing: setting up msg_gen for rt_none (graph generation)\n" << std::flush;
#endif
        typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));

#ifdef PRINT_DEBUG
	std::cerr << "europar: self stabilizing: done setting up msg_gen for rt_none (graph generation)\n" << std::flush;
#endif	    

	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
#ifdef PRINT_DEBUG
	std::cerr << "europar: self stabilizing: called distribute for rt_none (graph generation)\n" << std::flush;
#endif
      }
	
      if (routing == rt_hypercube) {

        typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
	  
      if (routing == rt_rook) {

        typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	    
	distribute
	  (Graph500Iter(scale, e_start, a, b),
	   Graph500Iter(scale, e_start + e_count, a, b),
	   flip_pair(), std::back_inserter(edges),
	   distrib, trans, msg_gen);
      }
    }

#ifdef BGP_REPORT_MEMORY
    long long memSize, usedMem; 
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
    getMemSize(&memSize);
    getUsedMem(&usedMem);
    std::cout << trans.rank() << ": Edge distribution complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

#ifdef PRINT_DEBUG
    std::cerr << "europar: self stabilizing: about to construct the graph\n" << std::flush;
#endif

    Digraph g(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), 
	      make_generator_iterator(gen, uniform_int<weight_type>(1, C)), n, trans, distrib);

#ifdef PRINT_DEBUG
    std::cerr << "europar: self stabilizing: graph has been constructed\n" << std::flush;
#endif

    // Clear edge array above
    edges.clear();

    // Generate sources
    boost::uniform_int<uint64_t> rand_vertex(0, n-1);
      
    { amplusplus::scoped_epoch epoch(trans); }

    time_type gen_end = get_time();

    if (trans.rank() == 0)
      std::cout << "Graph generation took " << print_time(gen_end - gen_start) << "s\n";

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
    if (trans.rank() == 0) {
#endif 
    getMemSize(&memSize);
    getUsedMem(&usedMem);
    std::cout << trans.rank() << ": Graph construction complete -- " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
    }
#endif 
#endif

#ifdef PRINT_ET
    if (trans.rank() == 0)
      std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

    // Weight map
    typedef property_map<Digraph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = get(&WeightedEdge::weight, g);
    
    // Distance map
    std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
    typedef iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap>  DistanceMap;
    DistanceMap distance(distanceS.begin(), get(vertex_index, g));

    // Predecessor map
    std::vector<Vertex> predecessorS(num_vertices(g), graph_traits<Digraph>::null_vertex());
    typedef iterator_property_map<std::vector<Vertex>::iterator, VertexIndexMap>  PredecessorMap;
    PredecessorMap predecessor(predecessorS.begin(), get(vertex_index, g));

#ifdef PRINT_DEBUG
    if(trans.rank() == 0)
      std::cout << "About to enter the thread loop" << std::endl;
#endif

    for (unsigned int level = levels_min; level <= levels_max; level += levels_step) {
    for (weight_type delta = delta_min; delta <= delta_max; delta += delta_step) {
    for (unsigned int i = 0, num_threads = thread_num_vals[i]; i < thread_num_vals.size(); ++i, num_threads = thread_num_vals[i]) {
	time_type total_self_stabilizing_time = 0;
	std::pair<unsigned long long, unsigned long long> cache_stats(0,0);

	if (no_reductions) {
	  //
	  // self-stabilizing shortest paths no reductions
	  //

	  total_self_stabilizing_time = 0;
	  cache_stats = std::make_pair(0, 0);

	  if (routing == rt_none) {

	    typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	    MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_size)));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
	        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

#ifdef PRINT_DEBUG
	      if(trans.rank() == 0)
		std::cout << "We are here" << std::endl;
#endif
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

#ifdef PRINT_DEBUG
	      std::cout << "About to start self-stabilizing shortest paths..." << std::endl;
#endif

	      time_type ss_time = run_self_stabilizing<Digraph, WeightMap, DistanceMap, PredecessorMap, MessageGenerator>
	        (trans, barrier_trans, g, weight, distance, predecessor, current_source, num_threads, n, verify, stats, msg_gen, ordering_size, level, delta);

	      if (ss_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

	      total_self_stabilizing_time += ss_time;
	    }
	  }

	  if (routing == rt_hypercube) {

	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
	        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ss_time = run_self_stabilizing<Digraph, WeightMap, DistanceMap, PredecessorMap, MessageGenerator>
	        (trans, barrier_trans, g, weight, distance, predecessor, current_source, num_threads, n, verify, stats, msg_gen, ordering_size, level, delta);

	      if (ss_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

	      total_self_stabilizing_time += ss_time;
	    }
	  }

	  if (routing == rt_rook) {

	    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
	        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ss_time = run_self_stabilizing<Digraph, WeightMap, DistanceMap, PredecessorMap, MessageGenerator>
	        (trans, barrier_trans, g, weight, distance, predecessor, current_source, num_threads, n, verify, stats, msg_gen, ordering_size, level, delta);

	      if (ss_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

	      total_self_stabilizing_time += ss_time;
	    }
	  }

	  if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	    unsigned long long buffers{}, messages{};
	    for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	      buffers += cumulative_flushes[i];
	      messages += cumulative_flushes[i] * i;
	    }
#endif
	    std::cout << "Total Self-Stabilizing (no reductions) time for scale = " << scale << ", for " << num_sources << " sources = " 
		      << print_time(total_self_stabilizing_time) << " (" 
		      << print_time(total_self_stabilizing_time / num_sources) << " per source)"
		      << ", " << num_threads << " threads, " << "ordering_size is " << ordering_size << ", coalescing_size = " << coalescing_size << ", delta = " << delta << ", levels = " << level
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		      << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		      << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		      << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		      << (messages/buffers) << " messages per buffer on average."
#endif
		      << std::endl;
	  }

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 
          getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif
	}

	if (per_thread_reductions) {
	//
	// Self-Stabilizing with reductions
	//

#ifdef PRINT_DEBUG
	  std::cerr << "europar: Entered per thread reductions\n" << std::flush;
#endif

	  total_self_stabilizing_time = 0;
	  cache_stats = std::make_pair(0, 0);

	  if (routing == rt_none) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::no_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
	        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ss_time = run_self_stabilizing<Digraph, WeightMap, DistanceMap, PredecessorMap, MessageGenerator>
	        (trans, barrier_trans, g, weight, distance, predecessor, current_source, num_threads, n, verify, stats, msg_gen, ordering_size, level, delta);

	      if (ss_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

	      total_self_stabilizing_time += ss_time;
	    }
	  }

	  if (routing == rt_hypercube) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	    
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
	        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif

	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ss_time = run_self_stabilizing<Digraph, WeightMap, DistanceMap, PredecessorMap, MessageGenerator>
	        (trans, barrier_trans, g, weight, distance, predecessor, current_source, num_threads, n, verify, stats, msg_gen, ordering_size, level, delta);

	      if (ss_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

	      total_self_stabilizing_time += ss_time;
	    }
	  }

	  if (routing == rt_rook) {

	    typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_size), reduction_cache_size, amplusplus::rook_routing(trans.rank(), trans.size()));
	
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
      clear_cumulative_work_stats();
#endif
	    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
	      if (trans.rank() == 0)
	        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
	  
	      Vertex current_source = vertex(rand_vertex(synch_gen), g);

	      time_type ss_time = run_self_stabilizing<Digraph, WeightMap, DistanceMap, PredecessorMap, MessageGenerator>
	        (trans, barrier_trans, g, weight, distance, predecessor, current_source, num_threads, n, verify, stats, msg_gen, ordering_size, level, delta);

	      if (ss_time == -1.) { // Not enough vertices visited
	        --source_i; continue;
	      } 

	      total_self_stabilizing_time += ss_time;
	    }
	  }
	
	  if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
	    unsigned long long buffers{}, messages{};
	    for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	      buffers += cumulative_flushes[i];
	      messages += cumulative_flushes[i] * i;
	    }
#endif
            std::cout << "Total Self-Stabilizing (w/ per-thread reductions) time for scale=" << scale << ", for " << num_sources << " sources = " 
		      << print_time(total_self_stabilizing_time) << " (" 
		      << print_time(total_self_stabilizing_time / num_sources) << " per source)"
		      << ", " << num_threads << " threads, " << "ordering_size is " << ordering_size << ", coalescing_size = " << coalescing_size << ", delta = " << delta << ", levels = " << level 
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
		      << "\nTests = " << cache_stats.first << ", hits = " << cache_stats.second
		      << ", hit rate = " << (double(cache_stats.second) / double(cache_stats.first) * 100.) << "%"
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS 
		      << ", messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), " 
		      << " full buffers = " << cumulative_full << " (" << (cumulative_full / num_sources) << " per source)" 
		      << " partial flushes = " << buffers << " (" << (buffers / num_sources) << " per source) with "
		      << (messages/buffers) << " messages per buffer on average."
#endif
		      << std::endl;
	  }

#ifdef BGP_REPORT_MEMORY
#if (BGP_REPORT_MEMORY == 1)
	  if (trans.rank() == 0) {
#endif 
          getMemSize(&memSize);
	  getUsedMem(&usedMem);
	  std::cout << trans.rank() << ": " << memSize << " bytes available (" << usedMem << " bytes in use)\n";
#if (BGP_REPORT_MEMORY == 1)
	  }
#endif 
#endif
	} // end if (per_thread_reductions)

    } // end for thread_num loop
    } // end levels loop
    } // end delta loop
  } // end all self-stabilizing shortest paths algorithms

  return 0;
}
