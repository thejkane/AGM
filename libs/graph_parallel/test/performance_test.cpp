// Copyright (C) 2014 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine
//           Marcin Zalewski
//           Thejaka Kanewala

#define LIB_CDS 1

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <inttypes.h>
#include <cstdlib>
#include <math.h>

//#define DISABLE_SELF_SEND_CHECK
//#define AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
//#define PBGL2_PRINT_WORK_STATS
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
//#define NUMBER_COMPONENTS
// #define PBGL_TIMING // turn on local component counting, sequential timing,  and such in CC
// #define PRINT_ET
#define BARRIER_TRANS
// #define CLONE
#define CALCULATE_GTEPS 1
#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>
#include <am++/counter_coalesced_message_type.hpp>
#include <boost/graph/distributed/time_calc.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/use_mpi.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/recursive_rmat_generator.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/graph500_generator.hpp>
#include <boost/graph/permute_graph.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths_node.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths_thread.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths_numa.hpp>
#include <boost/graph/distributed/mis.hpp>
#include <boost/graph/distributed/mis_delta.hpp>
#include <boost/graph/distributed/luby_mis.hpp>
#include <boost/graph/distributed/kla_sssp_buffer.hpp>
#include <boost/graph/distributed/kla_sssp_numa.hpp>
#include <boost/graph/distributed/kla_sssp_node.hpp>
#include <boost/graph/distributed/kla_sssp_thread.hpp>
#include <boost/graph/distributed/distributed_control.hpp>
#include <boost/graph/distributed/distributed_control_node.hpp>
#include <boost/graph/distributed/distributed_control_pheet.hpp>
#include <boost/graph/distributed/priority_q_defs.hpp>
#include <boost/graph/distributed/distributed_control_chaotic.hpp>
#include <boost/graph/distributed/cc_dc.hpp>
#include <boost/graph/distributed/cc_chaotic.hpp>
#include <boost/graph/distributed/delta_stepping_cc.hpp>
#include <boost/graph/distributed/cc_diff_dc.hpp>
#include <boost/graph/distributed/level_sync_cc.hpp>
#include <boost/graph/distributed/connected_components.hpp>
#include <boost/graph/distributed/sv_cc.hpp>
#include <boost/graph/distributed/page_rank.hpp>
#include <boost/graph/distributed/graph_reader.hpp>
#include <boost/graph/relax.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>

#include <boost/parallel/append_buffer.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <parallel/algorithm>
#include <algorithm> // for std::min, std::max
#include <functional>

#include <cds/init.h>       // for cds::Initialize and cds::Terminate
#include <cds/gc/hp.h>      // for cds::HP (Hazard Pointer) SMR

//#include "statistics.h"
#define HAS_NOT_BEEN_SEEN

#ifdef DISABLE_SELF_SEND_CHECK
#warning "==================^^^^^Self send check is disabled !^^^^^======================="
#endif

int _RANK;

using namespace boost;
using boost::parallel::all_reduce;
using boost::parallel::maximum;
using std::max;

#include <time.h>
#include <sys/time.h>

/*typedef double time_type;

inline time_type get_time()
{
  return MPI_Wtime();
#if 0
  timeval tp;
  gettimeofday(&tp, 0);
  return tp.tv_sec + tp.tv_usec / 1000000.0;
#endif
}*/

std::string print_time(time_type t)
{
  std::ostringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
  return out.str();
}

struct flip_pair {
  template <typename T>
  T operator()(const T& x) const {
    // flipping edges by keeping same weight
    return std::make_pair(x.second.first, std::make_pair(x.first, x.second.second));
  }
};


#define PRINT_GRAPH500_STATS(lbl, israte)				   \
  do {								   \
    printf ("Min(%s): %20.17e\n", lbl, stats[0]);		   \
    printf ("First Quartile(%s): %20.17e\n", lbl, stats[1]);	   \
    printf ("Median(%s): %20.17e\n", lbl, stats[2]);		   \
    printf ("Third Quartile(%s): %20.17e\n", lbl, stats[3]);	   \
    printf ("Max(%s): %20.17e\n", lbl, stats[4]);		   \
    if (!israte) {						   \
      printf ("Mean(%s): %20.17e\n", lbl, stats[5]);		   \
      printf ("Stddev(%s): %20.17e\n", lbl, stats[6]);		   \
      } else {							   \
      printf ("Harmonic Mean(%s): %20.17e\n", lbl, stats[7]);	   \
      printf ("Harmonic Stddev(%s): %20.17e\n", lbl, stats[8]);	   \
    }								   \
  } while (0)



extern void statistics (double *out, double *data, size_t n);

static int
dcmp (const void *a, const void *b)
{
  const double da = *(const double*)a;
  const double db = *(const double*)b;
  if (da > db) return 1;
  if (db > da) return -1;
  if (da == db) return 0;
  fprintf (stderr, "No NaNs permitted in output.\n");
  abort ();
  return 0;
}

void statistics (double *out, double *data, size_t n)
{
  long double s, mean;
  double t;
  int k;

  /* Quartiles */
  qsort (data, n, sizeof (*data), dcmp);
  out[0] = data[0];
  t = (n+1) / 4.0;
  k = (int) t;
  if (t == k)
    out[1] = data[k];
  else
    out[1] = 3*(data[k]/4.0) + data[k+1]/4.0;
  t = (n+1) / 2.0;
  k = (int) t;
  if (t == k)
    out[2] = data[k];
  else
    out[2] = data[k]/2.0 + data[k+1]/2.0;
  t = 3*((n+1) / 4.0);
  k = (int) t;
  if (t == k)
    out[3] = data[k];
  else
    out[3] = data[k]/4.0 + 3*(data[k+1]/4.0);
  out[4] = data[n-1];

  s = data[n-1];
  for (k = n-1; k > 0; --k)
    s += data[k-1];
  mean = s/n;
  out[5] = mean;
  s = data[n-1] - mean;
  s *= s;
  for (k = n-1; k > 0; --k) {
    long double tmp = data[k-1] - mean;
    s += tmp * tmp;
  }
  out[6] = sqrt (s/(n-1));

  s = (data[0]? 1.0L/data[0] : 0);
  for (k = 1; k < n; ++k)
    s += (data[k]? 1.0L/data[k] : 0);
  out[7] = n/s;
  mean = s/n;

  /*
    Nilan Norris, The Standard Errors of the Geometric and Harmonic
    Means and Their Application to Index Numbers, 1940.
    http://www.jstor.org/stable/2235723
  */
  s = (data[0]? 1.0L/data[0] : 0) - mean;
  s *= s;
  for (k = 1; k < n; ++k) {
    long double tmp = (data[k]? 1.0L/data[k] : 0) - mean;
    s += tmp * tmp;
  }
  s = (sqrt (s)/(n-1)) * out[7] * out[7];
  out[8] = s;
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
std::unique_ptr<atomic_flush_type[]> flushes;
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
typedef std::tuple<unsigned long long, unsigned long long, unsigned long long, unsigned long long, unsigned long long, unsigned long long> all_work_stats_t;

// Distributed control work stats
all_work_stats_t dc_stats, ds_stats;

void clear_cumulative_work_stats() {
  dc_stats = all_work_stats_t{ 0ul, 0ul, 0ul, 0ul, 0ul, 0ul };
  ds_stats = all_work_stats_t{ 0ul, 0ul, 0ul, 0ul, 0ul, 0ul };
}

void print_q_stats(unsigned long long max_q_sz,
		   unsigned long long avg_max_q_size = 0) {
  unsigned long long max, min;
  MPI_Allreduce(&max_q_sz, &max, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_q_sz, &min, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);

  if (_RANK == 0)
    std::cout << "Maximum max q size : " << 
      max << " minimum max q size : " << min << std::endl;

  if (avg_max_q_size != 0) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    unsigned long long avg_max_sum;
    MPI_Allreduce(&avg_max_q_size, &avg_max_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (_RANK == 0)
      std::cout << "Average maximum queue size : " << (avg_max_sum / world_size) << std::endl; 
  }
}



void print_and_accumulate_work_stats(work_stats_t stats, all_work_stats_t &cummulative_stats, size_t edges_traverse) {
  work_stats_t temp_stats;
  MPI_Allreduce(&std::get<0>(stats), &std::get<0>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&std::get<1>(stats), &std::get<1>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&std::get<2>(stats), &std::get<2>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&std::get<3>(stats), &std::get<3>(temp_stats), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(_RANK == 0) std::cout << "Useful work: " << std::get<0>(temp_stats) << ", invalidated work: " << std::get<1>(temp_stats) << ", useless work: " << std::get<2>(temp_stats) << ", rejected work: " << std::get<3>(temp_stats) << ", number of edges traversed : " << edges_traverse <<  ", level of ordering: " << (std::get<0>(temp_stats) + std::get<2>(temp_stats) - (2*edges_traverse)-1) << ", extra edge relaxes: " << (std::get<1>(temp_stats)-std::get<3>(temp_stats)) << std::endl;
  std::get<0>(cummulative_stats) += std::get<0>(temp_stats);
  std::get<1>(cummulative_stats) += std::get<1>(temp_stats);
  std::get<2>(cummulative_stats) += std::get<2>(temp_stats);
  std::get<3>(cummulative_stats) += std::get<3>(temp_stats);
  std::get<4>(cummulative_stats) += (std::get<0>(temp_stats) + std::get<2>(temp_stats) - (2*edges_traverse) - 1); // ordering
  std::get<5>(cummulative_stats) += (std::get<1>(temp_stats)-std::get<3>(temp_stats)); // extra relaxes
  
}
// TODO: Unify these two macros?  We probably don't need macros at all.  This is just a leftover from the previous test harness.
#define PBGL2_DC_PRINT \
  << "\nUseful work: " << (std::get<0>(dc_stats) / num_sources) << " (per source), invalidated work: " << (std::get<1>(dc_stats) / num_sources) << " (per source), useless work: " << (std::get<2>(dc_stats)/num_sources) << " (per source), rejected work: " << (std::get<3>(dc_stats)/num_sources) << " (per source). Rejected/useful ratio: " << ((double)std::get<3>(dc_stats) / (double)std::get<0>(dc_stats)) << ". Invalidated/useful ratio: " << ((double)std::get<1>(dc_stats) / (double)std::get<0>(dc_stats)) << ", level of ordering: " << std::get<4>(dc_stats) << ", extra edge relaxes: " << std::get<5>(dc_stats)
#define PBGL2_DS_PRINT \
  << "\nUseful work: " << (std::get<0>(ds_stats) / num_sources) << " (per source), invalidated work: " << (std::get<1>(ds_stats) / num_sources) << " (per source), useless work: " << (std::get<2>(ds_stats)/num_sources) << " (per source), rejected work: " << (std::get<3>(ds_stats)/num_sources) << " (per source). Rejected/useful ratio: " << ((double)std::get<3>(ds_stats) / (double)std::get<0>(ds_stats)) << ". Invalidated/useful ratio: " << ((double)std::get<1>(ds_stats) / (double)std::get<0>(ds_stats)) << ", level of ordering: " << (std::get<4>(ds_stats)/num_sources) << ", extra edge relaxes: " << (std::get<5>(ds_stats)/num_sources)
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


template <typename Graph, typename DistanceMap, typename WeightMap>
size_t get_gteps(amplusplus::transport& trans,  Graph& g, DistanceMap& distance, const WeightMap& weight) {
  distance.set_consistency_model(boost::parallel::cm_forward);
  distance.set_max_ghost_cells(0);

  size_t gteps = 0;
  {
    amplusplus::scoped_epoch epoch(g.transport());
    size_t dist = std::numeric_limits<size_t>::max();
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      dist = get(distance, v);
      if(dist < std::numeric_limits<size_t>::max()) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  gteps+=1;
	  //get(distance, target(e, g));
	}
      }
    }
  }
  gteps = gteps / 2; // we are considering an undirected graph
  return gteps;
}

template <typename Graph, typename DistanceMap, typename WeightMap, typename Vertex>
size_t count_edges(amplusplus::transport& trans,  Graph& g, DistanceMap& distance, const WeightMap& weight, const Vertex source) {
  distance.set_consistency_model(boost::parallel::cm_forward);
  distance.set_max_ghost_cells(0);

  size_t gteps = 0;
  {
    amplusplus::scoped_epoch epoch(g.transport());
    size_t dist = std::numeric_limits<size_t>::max();
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      dist = get(distance, v);
      if(dist < std::numeric_limits<size_t>::max()) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  // check whether this is a self loop; if self loop dont count it
	  if (target(e,g) == v)
	    continue;

	  gteps+=1;
	  //get(distance, target(e, g));
	}
      }

    }
  }
  return gteps;
}

template<typename MessageType>
class threaded_distributor {
public:
  threaded_distributor(MessageType& mtype) : msg_type(mtype) {}
  void operator()() {
  }
  template <typename InIter, typename Flip>
  void operator()(int tid, InIter begin, InIter end, const Flip& flip,
                  amplusplus::transport& trans) {
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      concurrent_distribute(begin, end, flip, trans, msg_type);
    }
  }

private:
  MessageType& msg_type;
};

template <typename Edges>
class threaded_merger {
public:
  threaded_merger(Edges& all) : all_edges(all) {}
  void operator()() {
  }
  template <typename EdgesIter>
  void operator()(EdgesIter pos, EdgesIter begin, EdgesIter end) {
    std::copy(begin, end, pos);
  }

private:
  Edges& all_edges;

};

template <typename Graph, typename ComponentMap, typename vertices_sz_t>
void calculate_cc_stats(Graph& g, ComponentMap& components, vertices_sz_t n, bool allstats) {
  typedef unsigned int component_t;
  typedef std::map<component_t, component_t> component_stats_t;
  component_stats_t component_stats;

#ifdef PRINT_DEBUG
  std::cout << "Calculating local components ..." << std::endl;
#endif

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    component_t vcomp = (component_t)get(components, v);
    component_stats_t::iterator 
      ite = component_stats.find(vcomp);
    if (ite == component_stats.end()) {
      component_stats.insert(std::make_pair(vcomp, 1));
    } else {
      ++((*ite).second);
    }
  }

  // end calculating local components.

  // exchange all components
  component_t localcomps = component_stats.size();
#ifdef PRINT_DEBUG
  std::cout << "Rank : " << _RANK << ", found local components : " << localcomps << std::endl;
#endif

  // allstats false
  if (!allstats) {
    // sort by number of vertices
    std::vector<std::pair<int, int>> pairs;
    for (auto itr = component_stats.begin(); itr != component_stats.end(); ++itr)
      pairs.push_back(*itr);

    component_stats.clear();

    std::sort(pairs.begin(), pairs.end(), [=](std::pair<int, int>& a, std::pair<int, int>& b) {
	return a.second > b.second;
      });

    if (_RANK == 0) {
      std::cout << "Printing component statistics ..." << std::endl;
      std::cout << "Printing only five largest components. To print all use --allstats" << std::endl;
      std::cout << "Total number of components in root : " << pairs.size() << std::endl;
    }

    int i = 0;
    vertices_sz_t rn = 0;
    for (; i < 5; ++i) {
      int compid = pairs[i].first;
      int vtcs = pairs[i].second;

      //check whether component id is same across all nodes
      int retcomp;
      MPI_Barrier(MPI_COMM_WORLD);
      std::cout << "r : " << _RANK << ", c=" << compid << ", i=" << i << std::endl;
      MPI_Allreduce(&compid, &retcomp, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      if (retcomp == compid) {
	// run a all reduce on all the vertices
	int allvs;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&vtcs, &allvs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (_RANK == 0) {
	  std::cout << "component : " << compid << " vertices : " << allvs
		    << " percentage : " << ((((double)allvs)/(double)n) * 100)
		    << std::endl;
	}

      } else {
	std::cout << "[ComponentStats] A mismatch found" << std::endl;
      }
    }
  } else { // allstats true

    component_t allcomps;
    // calculate sum of localcomps
    MPI_Allreduce(&localcomps, &allcomps, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    component_t expectedcomps = allcomps - localcomps;

    if (_RANK == 0)
      std::cout << "All the components : " << allcomps << ", expecting receive components : "
		<< expectedcomps << std::endl;

    typedef struct component_vertices {
      component_t c;
      component_t vs;

      component_vertices():
	c(0), vs(0) {}

      component_vertices(int pc, int pvs):
	c(pc), vs(pvs) {}
      component_vertices(const component_vertices &cvob):
	c(cvob.c), vs(cvob.vs) {}
    } cvs;


    const int nitems=2;
    int          blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_UNSIGNED, MPI_UNSIGNED};
    MPI_Datatype mpi_cv_type;
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(cvs, c);
    offsets[1] = offsetof(cvs, vs);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cv_type);
    MPI_Type_commit(&mpi_cv_type);

    std::vector<cvs> sendcomps;
    std::vector<cvs> recvcomps;

    if (_RANK != 0) {
      // go through local component stat maps and collect them to an array
      sendcomps.resize(localcomps);
      int j = 0;
      component_stats_t::iterator itesend = component_stats.begin();
      for (; itesend != component_stats.end(); ++itesend) {
	sendcomps[j].c = (*itesend).first;
	sendcomps[j].vs = (*itesend).second;
#ifdef PRINT_DEBUG
	std::cout << "Rank : " << _RANK << ", c=" << sendcomps[j].c << ", vs=" << sendcomps[j].vs << std::endl;
#endif
	++j;
      }

#ifdef PRINT_DEBUG
      std::cout << "Rank : " << _RANK << ", sending components : " << localcomps << std::endl;
#endif
      if (MPI_Send(&sendcomps[0], localcomps, mpi_cv_type, 0/*rank = 0*/, 13 /*tag*/,
		   MPI_COMM_WORLD) != 0) {
	std::cout << "MPI_Send : " << "Error exchanging components" << std::endl;
      }
    } 

    MPI_Barrier(MPI_COMM_WORLD);

    if (_RANK == 0) {
      recvcomps.resize(expectedcomps);
    
      if (expectedcomps != 0) {
#ifdef PRINT_DEBUG
	std::cout << "Rank : " << _RANK << ", expected components : " << expectedcomps << std::endl;
#endif
	int recvcount = 0;
	component_t rcvexpectedcomps = expectedcomps; 
	while (recvcount != expectedcomps) {
	  MPI_Status status;
	  if (MPI_Recv(&recvcomps[recvcount], rcvexpectedcomps, mpi_cv_type, MPI_ANY_SOURCE, 13,
		       MPI_COMM_WORLD, &status) != 0)
	    std::cout << "MPI_Recv : " << "Error receiving messages" << std::endl;
	
	  int instcount = 0;
	  if (MPI_Get_count(&status, mpi_cv_type, &instcount) != MPI_SUCCESS)
	    std::cout << "MPI_Get_count :" << "Error getting received count" << std::endl;

	  recvcount += instcount;
	  rcvexpectedcomps -= instcount;
	}
      }
    }

    sendcomps.clear();

#ifdef PRINT_DEBUG
    std::cout << "Rank : " << _RANK 
	      << ", Merging received components to component stats. Received components : "
	      << recvcomps.size() << std::endl;
#endif
    // update compnent stats
    typename std::vector<cvs>::iterator veciter = recvcomps.begin();
    for(; veciter != recvcomps.end(); ++veciter) {
      component_stats_t::iterator 
	ite = component_stats.find((*veciter).c);
      if (ite == component_stats.end()) {
#ifdef PRINT_DEBUG
	std::cout << "inserting c=" << (*veciter).c << ", vs=" << (*veciter).vs << std::endl; 
#endif
	component_stats.insert(std::make_pair((*veciter).c, (*veciter).vs));
      } else {
#ifdef PRINT_DEBUG
	std::cout << "adding c=" << (*veciter).c << ", vs=" << (*veciter).vs << std::endl;
#endif
	(*ite).second += (*veciter).vs;
      }
    }

    recvcomps.clear();


    if (_RANK == 0) {
      // sort by number of vertices
      std::vector<std::pair<int, int>> pairs;
      for (auto itr = component_stats.begin(); itr != component_stats.end(); ++itr)
	pairs.push_back(*itr);

      component_stats.clear();

#ifdef PRINT_DEBUG
      std::cout << "Sorting components by the number of vertices in their" << std::endl;
#endif
      std::sort(pairs.begin(), pairs.end(), [=](std::pair<int, int>& a, std::pair<int, int>& b) {
	  return a.second > b.second;
	});

      std::cout << "Printing component statistics ..." << std::endl;

      std::cout << "Total number of components : " << pairs.size() << std::endl;
      int i = 0;
      vertices_sz_t rn = 0;
      std::vector<std::pair<int, int>>::iterator piter = pairs.begin();
      for (; piter != pairs.end(); ++piter) {
	std::cout << "component : " << (*piter).first << " vertices : " << (*piter).second
		  << " percentage : " << ((((double)((*piter).second))/(double)n) * 100)
		  << std::endl;
	rn += (vertices_sz_t)(*piter).second;
      }


      std::cout << "rn = " << rn << ", n = " << n << std::endl;
      assert(rn == n);


      std::cout << "Done with CC stats." << std::endl;
    }
  }
}


template <typename Graph>
bool run_grah_stats(Graph& g) {
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  OwnerMap owner = get(vertex_owner, g);
  int remote_edges = 0;
  int local_edges = 0;
  BGL_FORALL_EDGES_T(e, g, Graph) {
    if (get(owner, source(e, g)) != _RANK)
      std::cout << "sv : " << source(e, g) << " svown : " 
		<< get(owner, source(e, g))
		<< " tv : " << target(e, g)
		<< " tvown : " << get(owner, target(e, g))
		<< " rank : " << _RANK << std::endl;

    //    assert(get(owner, source(e, g)) == _RANK);
    if (get(owner, target(e, g)) != _RANK)
      ++remote_edges;
    else
      ++local_edges;
  }

  int alllocale = num_edges(g);
  std::cout << "Rank : " << _RANK << " remote e : " << remote_edges
	    << " local e : " << local_edges 
	    << " remote per : " 
	    << ((double)remote_edges / (double)alllocale) * 100 
	    << "%" << std::endl;

}

template <typename Graph, typename ComponentMap>
bool verify_cc(Graph& g, ComponentMap& components) {
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  components.set_consistency_model(boost::parallel::cm_forward);
  components.set_max_ghost_cells(0);
	      
  {
    amplusplus::scoped_epoch epoch(g.transport());
      
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	get(components, source(e, g));
	get(components, target(e, g));
      }
    }
  }

  {	    
    amplusplus::scoped_epoch epoch(g.transport()); // at the moment get() sends a message
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_ADJ_T(v, u, g, Graph) {
	//	  std::cout << "verifying vertex v : " << v << std::endl;
	//#ifdef PRINT_DEBUG
	if (get(components, v) != get(components, u)) 
	  std::cout << "Component of " << v << " : " << get(components, v)
		    << " component of " << u << " : " << get(components, u)
		    << std::endl;

	assert(get(components, v) == get(components, u)); 
      }
    }
  }
    
  components.clear(); // Clear memory used by ghost cells
}

template <typename Graph, typename DistanceMap, typename WeightMap>
bool verify_sssp(amplusplus::transport& trans,  Graph& g, DistanceMap& distance, const WeightMap& weight) {
  distance.set_consistency_model(boost::parallel::cm_forward);
  distance.set_max_ghost_cells(0);

  if (trans.rank() == 0) std::cout<<"Verifying results......";

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
      if(get(distance, target(e, g)) > boost::closed_plus<weight_type>()(get(distance, v), get(weight, e))) std::abort();
#endif
    }
  }
      if (trans.rank() == 0) std::cout << "Verified." << std::endl;
      distance.clear(); // Clear memory used by ghost cells

      return true;
}

template <typename Graph, typename MISMap>
bool verify_mis(amplusplus::transport& trans,  Graph& g, MISMap& mis) {
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  OwnerMap owner(get(vertex_owner, g));

  mis.set_consistency_model(boost::parallel::cm_forward || boost::parallel::cm_backward);
  mis.set_max_ghost_cells(0);

  if (trans.rank() == 0) std::cout<<"Verifying MIS results......";

  {
    amplusplus::scoped_epoch epoch(g.transport());
	  
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	get(mis, target(e, g));
      }
    }
  }
	    
  bool result = true;
  int found_mis = 0;
#ifdef PRINT_DEBUG
  std::cout << "MIS = {";
#endif
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (get(mis, v) ==  MIS_FIX1) {// v in mis, none of the neigbours should be in mis
#ifdef PRINT_DEBUG
      if (v == 35) {
	  std::cout << "Printing all neighbors of " << v
		    << " -- [";

	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  std::cout << target(e, g) << ", ";
	}
	std::cout << "]" << std::endl;
      }

      std::cout << v << "-[";
#endif
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
#ifdef PRINT_DEBUG
	std::cout << target(e, g) << ", ";
#endif
	if (v == target(e, g))
	  continue;

	if (get(mis, target(e, g)) == MIS_FIX1) {
	  std::cout << "[FAIL] Vertex : " << v << " and its neigbour " << target(e, g)
		    << " are both in MIS." 
		    << " Vertex value : " << get(mis, v)
		    << " neighbor value : " << get(mis, target(e, g)) 
		    << std::endl;
	  std::cout << "Neighbors of " << v << "- {";
	  BGL_FORALL_ADJ_T(v, u1, g, Graph) {
	    std::cout << "(" << u1 << ", " << get(mis, u1) << ")";
	  }
	  std::cout << "}" << std::endl;

	  auto k = target(e, g);
	  if (get(owner, k) == g.transport().rank()) {
	    std::cout << "Neighbors of " << k << "- {";
	    BGL_FORALL_ADJ_T(k, u2, g, Graph) {
	      std::cout << "(" << u2 << ", " << get(mis, u2) << ")";
	    }
	    std::cout << "}" << std::endl;
	  }

	  result = false;
	} else {
#ifdef PRINT_DEBUG
	  if (get(mis, target(e, g)) != MIS_FIX0) {
	    std::cout << "vertex : " 
		      << target(e, g) 
		      << ", is in " 
		      << get(mis, target(e, g))
		      << std::endl;
	  }
#endif

	  // cannot be MIS_UNFIX
	  assert(get(mis, target(e, g)) == MIS_FIX0);
	  found_mis = 1;
	}
	}
#ifdef PRINT_DEBUG
      std::cout << "], ";
#endif
    } else {
      // if v is not in mis, at least one of its neighbors must
      // be in mis
      if (get(mis, v) != MIS_FIX0) {
	std::cout << "[ERROR] Vertex - " << v << " is in " << get(mis, v)
		  << std::endl;
      }

      if (get(mis, v) != MIS_FIX0) {
	std::cout << "Error : " << v << " is not in MIS_FIX0. But in " << get(mis, v)
		  << std::endl;
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  auto k = target(e, g);
	  std::cout << "neigbor : " << k << " mis : " << get(mis, k)
		    << std::endl;
	}
      }
      assert(get(mis, v) == MIS_FIX0);

      bool inmis = false;
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	if (v == target(e, g))
	  continue;

	if (get(mis, target(e, g)) == MIS_FIX1) {
	  inmis = true;
	  break;
	}
      }
      
      if (!inmis) {
	std::cout << "Vertex : " << v << " and none of its neighbors is in MIS" << std::endl;
	assert(false);
      }
    }
  }
#ifdef PRINT_DEBUG
  std::cout << "}" << std::endl;
#endif

  // Run a reduction on found mis.
  // We cannot expect every rank will have found_mis true.
  // E.g:- rank0 - {0,1}, rank1 - {281474976710656, 281474976710657}
  // If every vertex is connect to each other, a vertex (Lets say 0) will
  // be in one rank. In the second rank all vertices are out of mis.
  int or_found_mis = 0;
  MPI_Allreduce(&found_mis, &or_found_mis, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (or_found_mis == 0) {
    std::cout << "Did not find any MIS" <<std::endl;
    result = false;
  }

  return result;
}

//=================================================//
// vertex id distributions
enum id_distribution_t {
  vertical,
  horizontal
};
//=================================================//


template <typename PriorityQueueGenerator = boost::graph::distributed::thread_priority_queue_gen,
	  typename Graph, typename MISMap,
	  typename MessageGenerator>
time_type
run_fix_mis(amplusplus::transport& trans, 
		  amplusplus::transport& barrier_trans, 
		  Graph& g,  
		  MISMap& mis, 
		  int num_threads, 
		  typename graph_traits<Graph>::vertices_size_type n, 
		  bool verify, 
		  MessageGenerator msg_gen, 
		  int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;


  if (trans.rank() == 0)
    std::cout << "Initializing mis map ..." << std::endl;

  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(mis, v, MIS_UNFIX); }

  trans.set_nthreads(num_threads);

  if (trans.rank() == 0)
    std::cout << "Creating algorithm instance ..." << std::endl;

  boost::graph::distributed::maximal_independent_set<Graph, MISMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, mis, trans, flushFreq, msg_gen);
    
  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
  // Many threads now
  trans.set_nthreads(num_threads);

  if (trans.rank() == 0)
    std::cout << "Invoking algorithm ..." << std::endl;

  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
	  
  D(0);
    
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
	  
  time_type end = get_time();

  if (trans.rank() == 0)
    std::cout << "Algorithm done ..." << std::endl;

  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

  //#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
  //#endif


  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    if (trans.rank()==0)
      std::cout << "Verifying mis ..." << std::endl;

    if (!verify_mis(trans, g, mis)) {
      std::cout << "MIS Verification Failed" << std::endl;
      assert(false);
      return 0;
    }
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(mis, v) != MIS_UNFIX)
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< std::endl;

  //if (total < 100) return -1.;
	  
  return end - start;
}

template <typename Graph, typename MISMap,
	  typename MessageGenerator>
time_type
run_fix_mis_bucket(amplusplus::transport& trans, 
		  amplusplus::transport& barrier_trans, 
		  Graph& g,  
		  MISMap& mis, 
		  int num_threads, 
		  typename graph_traits<Graph>::vertices_size_type n, 
		  bool verify, 
		  MessageGenerator msg_gen, 
		  int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;


  if (trans.rank() == 0)
    std::cout << "Initializing mis-delta map ..." << std::endl;

  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(mis, v, MIS_UNFIX); }

  trans.set_nthreads(num_threads);

  if (trans.rank() == 0)
    std::cout << "Creating algorithm instance ..." << std::endl;

  boost::graph::distributed::maximal_independent_set_delta<Graph, MISMap,
    append_buffer<Vertex, 10u>, MessageGenerator> 
    D(g, mis, trans, flushFreq, msg_gen);
    
  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
  // Many threads now
  trans.set_nthreads(num_threads);

  if (trans.rank() == 0)
    std::cout << "Invoking algorithm ..." << std::endl;

  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
	  
  D(0);
    
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
	  
  time_type end = get_time();

  if (trans.rank() == 0)
    std::cout << "Algorithm done ..." << std::endl;

  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

  //#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
  //#endif


  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    if (trans.rank()==0)
      std::cout << "Verifying delta mis ..." << std::endl;

    if (!verify_mis(trans, g, mis)) {
      std::cout << "Delta-MIS Verification Failed" << std::endl;
      assert(false);
      return 0;
    }
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(mis, v) != MIS_UNFIX)
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< std::endl;

  //if (total < 100) return -1.;
	  
  return end - start;
}


template <typename SelectGenerator, typename Graph, typename MISMap, typename MessageGenerator>
time_type
run_luby_maximal_is(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g,  
	     MISMap& mis, 
	     int num_threads, 
	     typename graph_traits<Graph>::vertices_size_type n, 
	     bool verify, 
	     MessageGenerator msg_gen, 
	     int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef iterator_property_map<typename std::vector<random_t>::iterator, VertexIndexMap>  RandomMap;

  // create a property map
  std::vector<random_t> pivec(num_vertices(g), 0);
  RandomMap rmap(pivec.begin(), get(vertex_index, g)); // TODO remove copying

  if (trans.rank() == 0)
    std::cout << "Initializing mis map ..." << std::endl;

  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(mis, v, MIS_UNFIX); }

  trans.set_nthreads(num_threads);

  if (trans.rank() == 0)
    std::cout << "Creating algorithm instance ..." << std::endl;

  boost::graph::distributed::luby_mis<Graph, MISMap, RandomMap, 
    //    boost::graph::distributed::select_a_functor_gen,
    SelectGenerator,
    MessageGenerator>
    D(g, mis, rmap, n, trans, flushFreq, msg_gen);
    
  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
  // Many threads now
  trans.set_nthreads(num_threads);

  if (trans.rank() == 0)
    std::cout << "Invoking algorithm ..." << std::endl;

  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
	  
  D(0);
    
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
	  
  time_type end = get_time();

  if (trans.rank() == 0)
    std::cout << "Algorithm done ..." << std::endl;

  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

  //#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
  //#endif


  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    if (trans.rank()==0)
      std::cout << "Verifying mis ..." << std::endl;

    if (!verify_mis(trans, g, mis)) {
      std::cout << "MIS Verification Failed" << std::endl;
      assert(false);
      return 0;
    }
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(mis, v) != MIS_UNFIX)
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< std::endl;

  //if (total < 100) return -1.;
	  
  return end - start;
}


template <typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_kla_sssp(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g, const WeightMap& weight, 
	     DistanceMap& distance, 
	     typename graph_traits<Graph>::vertex_descriptor current_source, 
	     int num_threads, 
	     typename graph_traits<Graph>::vertices_size_type n, 
	     bool verify, bool level_sync, 
	     MessageGenerator msg_gen, size_t k_level) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename property_traits<WeightMap>::value_type Dist;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);
    
  boost::graph::distributed::kla_shortest_paths_buffer<Graph, DistanceMap, WeightMap,
    append_buffer<std::pair<Vertex, std::pair<Dist, size_t> >, 10u>, MessageGenerator>
    D(g, distance, weight, trans, k_level, msg_gen);
    
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) {
    verify_sssp(trans, g, distance, weight); 
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}



template <typename PriorityQueueGenerator = boost::graph::distributed::numa_priority_queue_gen,
  typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_kla_sssp_numa(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g, const WeightMap& weight, 
	     DistanceMap& distance, 
	     typename graph_traits<Graph>::vertex_descriptor current_source, 
	     int num_threads, 
	     typename graph_traits<Graph>::vertices_size_type n, 
	     bool verify, bool level_sync, 
		  MessageGenerator msg_gen, size_t k_level,
		  int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename property_traits<WeightMap>::value_type Dist;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);

#ifdef LIB_CDS
  // LIBCDS Stuff goes here
  if (trans.rank() == 0)
    std::cout << "Initializing LIBCDS .... " << std::endl;

  cds::Initialize();
  // Initialize Hazard Pointer singleton
  cds::gc::HP hpGC;
  // Attach for the main thread
  cds::threading::Manager::attachThread();
#endif
    
  boost::graph::distributed::kla_shortest_paths_numa<Graph, DistanceMap, WeightMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, k_level, flushFreq, msg_gen);
    
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
  time_type start = D.get_start_time();

#ifdef LIB_CDS
  if (trans.rank() == 0)
    std::cout << "Termminating LIBCDS .... " << std::endl;

  cds::Terminate();
#endif

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) {
    verify_sssp(trans, g, distance, weight); 
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}

template <typename PriorityQueueGenerator = boost::graph::distributed::numa_priority_queue_gen,
  typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_kla_sssp_thread(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g, const WeightMap& weight, 
	     DistanceMap& distance, 
	     typename graph_traits<Graph>::vertex_descriptor current_source, 
	     int num_threads, 
	     typename graph_traits<Graph>::vertices_size_type n, 
	     bool verify, bool level_sync, 
		    MessageGenerator msg_gen, size_t k_level,
		    int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename property_traits<WeightMap>::value_type Dist;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);

  boost::graph::distributed::kla_shortest_paths_thread<Graph, DistanceMap, WeightMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, k_level,flushFreq, msg_gen);
    
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) {
    verify_sssp(trans, g, distance, weight); 
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}



template <typename PriorityQueueGenerator = boost::graph::distributed::numa_priority_queue_gen,
  typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_kla_sssp_node(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g, const WeightMap& weight, 
	     DistanceMap& distance, 
	     typename graph_traits<Graph>::vertex_descriptor current_source, 
	     int num_threads, 
	     typename graph_traits<Graph>::vertices_size_type n, 
	     bool verify, bool level_sync, 
		  MessageGenerator msg_gen, size_t k_level,
		  int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename property_traits<WeightMap>::value_type Dist;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);

#ifdef LIB_CDS
  // LIBCDS Stuff goes here
  if (trans.rank() == 0)
    std::cout << "Initializing LIBCDS .... " << std::endl;

  cds::Initialize();
  // Initialize Hazard Pointer singleton
  cds::gc::HP hpGC;
  // Attach for the main thread
  cds::threading::Manager::attachThread();
#endif
    
  boost::graph::distributed::kla_shortest_paths_node<Graph, DistanceMap, WeightMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, k_level, flushFreq, msg_gen);
    
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
        
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
  time_type start = D.get_start_time();

#ifdef LIB_CDS
  if (trans.rank() == 0)
    std::cout << "Termminating LIBCDS .... " << std::endl;

  cds::Terminate();
#endif

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) {
    verify_sssp(trans, g, distance, weight); 
  }

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}


template <typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_delta_stepping(amplusplus::transport& trans, amplusplus::transport& barrier_trans, Graph& g, const WeightMap& weight, DistanceMap& distance, typename graph_traits<Graph>::vertex_descriptor current_source, weight_type delta, int num_threads, typename graph_traits<Graph>::vertices_size_type n, bool verify, bool level_sync, MessageGenerator msg_gen) {
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
      D(g, distance, weight, trans, delta, msg_gen);
    trans.set_nthreads(1);

    D.set_source(current_source);
    
    if (level_sync)
      D.set_level_sync();

    { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
    time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    print_and_clear_epoch_times();
    print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
    const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
    print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
    size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
    work_stats_t work_stats = D.get_work_stats();
    print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

    // Back to one thread
    trans.set_nthreads(1);

    unsigned long long num_levels = D.get_num_levels();

    if (verify) verify_sssp(trans, g, distance, weight); 

    vertices_size_type visited = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (get(distance, v) < std::numeric_limits<weight_type>::max())
	++visited; 
    }
	  
    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
      r(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r(visited);
	  
    if (verify)
      if (trans.rank() == 0)
	std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		  << ",  " << num_levels << " levels required." << std::endl;

    if (total < 100) return -1.;
	  
    return end - start;
}


template <typename PriorityQueueGenerator = boost::graph::distributed::default_priority_queue_gen, 
	  typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_delta_stepping_node(amplusplus::transport& trans, 
		       amplusplus::transport& barrier_trans, 
		       Graph& g, const WeightMap& weight, 
		       DistanceMap& distance, 
		       typename graph_traits<Graph>::vertex_descriptor current_source, 
		       weight_type delta, int num_threads, typename graph_traits<Graph>::
			vertices_size_type n, bool verify, bool level_sync, MessageGenerator msg_gen,
			int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);

#ifdef LIB_CDS
  // LIBCDS Stuff goes here
  if (trans.rank() == 0)
    std::cout << "Initializing LIBCDS .... " << std::endl;

  cds::Initialize();
  // Initialize Hazard Pointer singleton
  cds::gc::HP hpGC;
  // Attach for the main thread
  cds::threading::Manager::attachThread();
#endif

  boost::graph::distributed::delta_stepping_shortest_paths_node<Graph, DistanceMap, WeightMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, delta, flushFreq, msg_gen);
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
  time_type start = D.get_start_time();

#ifdef LIB_CDS
  if (trans.rank() == 0)
    std::cout << "Termminating LIBCDS .... " << std::endl;

  cds::Terminate();
#endif

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) verify_sssp(trans, g, distance, weight); 

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}



template <typename PriorityQueueGenerator = boost::graph::distributed::default_priority_queue_gen, 
	  typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_delta_stepping_numa(amplusplus::transport& trans, 
		       amplusplus::transport& barrier_trans, 
		       Graph& g, const WeightMap& weight, 
		       DistanceMap& distance, 
		       typename graph_traits<Graph>::vertex_descriptor current_source, 
		       weight_type delta, int num_threads, typename graph_traits<Graph>::
			vertices_size_type n, bool verify, bool level_sync, MessageGenerator msg_gen,
			int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);

#ifdef LIB_CDS
  // LIBCDS Stuff goes here
  if (trans.rank() == 0)
    std::cout << "Initializing LIBCDS .... " << std::endl;

  cds::Initialize();
  // Initialize Hazard Pointer singleton
  cds::gc::HP hpGC;
  // Attach for the main thread
  cds::threading::Manager::attachThread();
#endif

  boost::graph::distributed::delta_stepping_shortest_paths_numa<Graph, DistanceMap, WeightMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, delta, flushFreq, msg_gen);
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
  time_type start = D.get_start_time();

#ifdef LIB_CDS
  if (trans.rank() == 0)
    std::cout << "Termminating LIBCDS .... " << std::endl;

  cds::Terminate();
#endif

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) verify_sssp(trans, g, distance, weight); 

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}



template <typename PriorityQueueGenerator = boost::graph::distributed::default_priority_queue_gen, 
	  typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_delta_stepping_thread(amplusplus::transport& trans, 
		       amplusplus::transport& barrier_trans, 
		       Graph& g, const WeightMap& weight, 
		       DistanceMap& distance, 
		       typename graph_traits<Graph>::vertex_descriptor current_source, 
		       weight_type delta, int num_threads, typename graph_traits<Graph>::
			  vertices_size_type n, bool verify, bool level_sync, MessageGenerator msg_gen,
			  int flushFreq) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  BGL_FORALL_VERTICES_T(v, g, Graph) 
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
	  
  trans.set_nthreads(num_threads);

  boost::graph::distributed::delta_stepping_shortest_paths_thread<Graph, DistanceMap, WeightMap,
    PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, delta, flushFreq, msg_gen);
  trans.set_nthreads(1);

  D.set_source(current_source);
    
  if (level_sync)
    D.set_level_sync();

  { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif
    
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
  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif

#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, ds_stats, edges_in_comp);
#endif

  // Back to one thread
  trans.set_nthreads(1);

  unsigned long long num_levels = D.get_num_levels();

  if (verify) verify_sssp(trans, g, distance, weight); 

  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited; 
  }
	  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
	  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start) 
		<< ",  " << num_levels << " levels required." << std::endl;

  if (total < 100) return -1.;
	  
  return end - start;
}


template <typename PriorityQueueGenerator = boost::graph::distributed::pheet_priority_queue_gen_128, typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_distributed_control_pheet(amplusplus::transport& trans, 
			     amplusplus::transport& barrier_trans, 
			     Graph& g, const WeightMap& weight, 
			     DistanceMap& distance, 
			     typename graph_traits<Graph>::vertex_descriptor current_source, 
			     int num_threads, typename graph_traits<Graph>::vertices_size_type n, 
			     bool verify, MessageGenerator msg_gen, 
			     MessageGenerator priority_msg_gen, 
			     unsigned int flushFreq, 
			     unsigned int eager_limit) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  // Initialize with infinite (max) distance
  BGL_FORALL_VERTICES_T(v, g, Graph)
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
  
  trans.set_nthreads(num_threads);

  boost::graph::distributed::distributed_control_pheet<Graph, DistanceMap, WeightMap, PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, msg_gen, priority_msg_gen, flushFreq, eager_limit);
  trans.set_nthreads(1);
  //trans.set_recvdepth(recvDepth);
  
  D.set_source(current_source);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
  
  D(0);
  
  // TODO: Add exception handling.  Excptions should be caught in every thread and checked in the main thread.
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();
  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif
  
#ifdef PBGL2_PRINT_WORK_STATS
  // calculate the number of edges
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, dc_stats, edges_in_comp);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);

#ifdef PBGL2_PRINT_WORK_STATS
  //get the queue size
  size_t local_q_size = D.get_max_q_size();
  print_q_stats(local_q_size);
#endif

  if (verify) verify_sssp(trans, g, distance, weight);
      
  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited;
  }
  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> >
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start)  << "\n" << std::endl;
  
  if (total < 100) return -1.;
  
  return end - start;
}


template <typename PriorityQueueGenerator = boost::graph::distributed::default_priority_queue_gen, typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_distributed_control_node(amplusplus::transport& trans, 
			     amplusplus::transport& barrier_trans, 
			     Graph& g, const WeightMap& weight, 
			     DistanceMap& distance, 
			     typename graph_traits<Graph>::vertex_descriptor current_source, 
			     int num_threads, typename graph_traits<Graph>::vertices_size_type n, 
			     bool verify, MessageGenerator msg_gen, 
			     MessageGenerator priority_msg_gen, 
			     unsigned int flushFreq, 
			     unsigned int eager_limit,
			     bool numa=false) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  // Initialize with infinite (max) distance
  BGL_FORALL_VERTICES_T(v, g, Graph)
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
  
  trans.set_nthreads(num_threads);

#ifdef LIB_CDS
  // LIBCDS Stuff goes here
  if (trans.rank() == 0)
    std::cout << "Initializing LIBCDS .... " << std::endl;

  cds::Initialize();
  // Initialize Hazard Pointer singleton
  cds::gc::HP hpGC;
  // Attach for the main thread
  cds::threading::Manager::attachThread();
#endif


  boost::graph::distributed::distributed_control_node<Graph, DistanceMap, WeightMap, PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, msg_gen, priority_msg_gen, flushFreq, eager_limit, numa);
  trans.set_nthreads(1);
  //trans.set_recvdepth(recvDepth);
  
  D.set_source(current_source);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
  
  D(0);
  
  // TODO: Add exception handling.  Excptions should be caught in every thread and checked in the main thread.
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();
  time_type start = D.get_start_time();

#ifdef LIB_CDS
  if (trans.rank() == 0)
    std::cout << "Termminating LIBCDS .... " << std::endl;

  cds::Terminate();
#endif
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif
  
#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, dc_stats, edges_in_comp);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);

#ifdef PBGL2_PRINT_WORK_STATS
  //get the queue size
  size_t local_q_size = D.get_max_q_size();
  print_q_stats(local_q_size);
#endif

  if (verify) verify_sssp(trans, g, distance, weight);
      
  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited;
  }
  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> >
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);
  
  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start)  << "\n" << std::endl;
  
  if (total < 100) return -1.;
  
  return end - start;
}

// Duplicating code; pretty bad ....
// TODO find the metaprogramming way to add cds initialization code
template <typename PriorityQueueGenerator = boost::graph::distributed::default_priority_queue_gen, typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_distributed_control(amplusplus::transport& trans, amplusplus::transport& barrier_trans, Graph& g, const WeightMap& weight, DistanceMap& distance, typename graph_traits<Graph>::vertex_descriptor current_source, int num_threads, typename graph_traits<Graph>::vertices_size_type n, bool verify, MessageGenerator msg_gen, MessageGenerator priority_msg_gen, unsigned int flushFreq, unsigned int eager_limit) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  // Initialize with infinite (max) distance
  BGL_FORALL_VERTICES_T(v, g, Graph)
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
  
  trans.set_nthreads(num_threads);

  boost::graph::distributed::distributed_control<Graph, DistanceMap, WeightMap, PriorityQueueGenerator, MessageGenerator>
    D(g, distance, weight, trans, msg_gen, priority_msg_gen, flushFreq, eager_limit);
  trans.set_nthreads(1);
  //trans.set_recvdepth(recvDepth);
  
  D.set_source(current_source);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  
  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
  
  D(0);
  
  // TODO: Add exception handling.  Excptions should be caught in every thread and checked in the main thread.
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();
  time_type start = D.get_start_time();

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
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  print_and_accumulate_work_stats(work_stats, dc_stats, edges_in_comp);

  // print q stats 
  print_q_stats(D.get_max_q_size(), D.get_avg_max_q_size());
#endif
  
  // Back to one thread
  trans.set_nthreads(1);
  
  if (verify) verify_sssp(trans, g, distance, weight);
      
  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited;
  }
  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> >
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);

  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start)  << "\n" << std::endl;
  
  if (total < 100) return -1.;
  
  return end - start;
}

// Duplicating code; pretty bad ....
// TODO find the metaprogramming way to add cds initialization code
template <typename Graph, typename WeightMap, typename DistanceMap, typename MessageGenerator>
time_type
run_distributed_control_chaotic(amplusplus::transport& trans, amplusplus::transport& barrier_trans, Graph& g, const WeightMap& weight, DistanceMap& distance, typename graph_traits<Graph>::vertex_descriptor current_source, int num_threads, typename graph_traits<Graph>::vertices_size_type n, bool verify, MessageGenerator msg_gen, unsigned int flushFreq, unsigned int eager_limit) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  
  // Initialize with infinite (max) distance
  BGL_FORALL_VERTICES_T(v, g, Graph)
    { put(distance, v, std::numeric_limits<weight_type>::max()); }
  
  trans.set_nthreads(num_threads);

  boost::graph::distributed::distributed_control_chaotic<Graph, DistanceMap, WeightMap, MessageGenerator>
    D(g, distance, weight, trans, msg_gen, flushFreq, eager_limit);
  trans.set_nthreads(1);
  //trans.set_recvdepth(recvDepth);
  
  D.set_source(current_source);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  
  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
  
  D(0);
  
  // TODO: Add exception handling.  Excptions should be caught in every thread and checked in the main thread.
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();
  time_type start = D.get_start_time();

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  print_and_accumulate_cache_stats(stats);
#endif
  
#ifdef PBGL2_PRINT_WORK_STATS
  size_t edges_in_comp =  count_edges(trans,  g, distance, weight, current_source);
  work_stats_t work_stats = D.get_work_stats();
  print_and_accumulate_work_stats(work_stats, dc_stats, edges_in_comp);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);
  
  if (verify) verify_sssp(trans, g, distance, weight);
      
  vertices_size_type visited = 0;
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    if (get(distance, v) < std::numeric_limits<weight_type>::max())
      ++visited;
  }
  
  boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> >
    r(trans, std::plus<vertices_size_type>());
  vertices_size_type total = r(visited);

  if (verify)
    if (trans.rank() == 0)
      std::cout << "Visited " << total << " vertices of " << n << " in " << print_time(end - start)  << "\n" << std::endl;
  
  if (total < 100) return -1.;
  
  return end - start;
}


template <typename Graph, typename MessageGenerator>
time_type
run_ps_sv_cc(amplusplus::transport& trans, 
	     Graph& g,
	     int num_threads, 
	     bool verify, 
	     bool level_sync, 
	     MessageGenerator msg_gen,
	     size_t vertices_per_lock) {

#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

  if (trans.rank() == 0) std::cout << "Initializing PS-SV-CC algorithm ...  " << std::endl;

  // Instantiate algorithms used later here so we don't time their constructions
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;

  // Parent, component, and lock maps for SV CC
  std::vector<Vertex> parentS(num_vertices(g), graph_traits<Graph>::null_vertex());
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap> ParentMap;
  ParentMap parent(parentS.begin(), get(vertex_index, g));

  std::vector<int> sv_componentS(num_vertices(g), std::numeric_limits<int>::max());
  typedef iterator_property_map<std::vector<int>::iterator, VertexIndexMap> ComponentMap;
  ComponentMap component(sv_componentS.begin(), get(vertex_index, g));

  typedef boost::parallel::lock_map<VertexIndexMap> LockMap;
  LockMap locks(get(vertex_index, g), num_vertices(g) / vertices_per_lock);

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(parent, v, v);
  }

  using boost::parallel::all_reduce;
  using boost::parallel::maximum;

  all_reduce<vertices_size_type, maximum<vertices_size_type> > 
    reduce_max(trans, maximum<vertices_size_type>());

  // TODO: If level_sync we should make the PS below a BFS

  boost::graph::distributed::parallel_search<Graph, ParentMap, MessageGenerator> 
    PS(g, parent, graph_traits<Graph>::null_vertex(), msg_gen);

  trans.set_nthreads(num_threads);
  boost::graph::distributed::connected_components<Graph, ParentMap, MessageGenerator>
    CC(trans, g, parent, locks, msg_gen);
  trans.set_nthreads(1);

  if (level_sync) 
    CC.set_level_sync();

  
  // Less than on vertices w/ special casing for null_vertex() so we start parallel 
  // search from the same vertex below
  typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  typedef typename property_map<Graph, vertex_local_t>::const_type LocalMap;
  typedef graph::distributed::cc_vertex_compare<OwnerMap, LocalMap> VertexLessThan;
  VertexLessThan vertex_lt(get(vertex_owner, g), get(vertex_local, g)); 

  //
  // Start timing loop
  //
  { amplusplus::scoped_epoch epoch(trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif

  if (trans.rank() == 0) std::cout << "Starting PS-SV-CC algorithm ...  " << std::endl;  	

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

  if (verify)
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

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif


  if (trans.rank() == 0) std::cout << "Done with PS-SV-CC algorithm in : " << (end-start) << std::endl;	
  // Back to one thread
  trans.set_nthreads(1);
  
#ifdef NUMBER_COMPONENTS
#error "Should be disabled for performance testing"
  int num_components = CC.template number_components<ComponentMap>(component);
  if (trans.rank() == 0)
    std::cout << num_components << " components found\n";
#endif

  if (verify) {
    if (trans.rank() == 0) std::cout << "Verifying PS-SV-CC algorithm ...  " << std::endl;

    parent.set_consistency_model(boost::parallel::cm_forward);
    parent.set_max_ghost_cells(0);
    
    {
      amplusplus::scoped_epoch epoch(g.transport());
      
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


template <typename Graph, typename MessageGenerator>
time_type
run_sv_c_c(amplusplus::transport& trans, Graph& g, 
	   int num_threads, 
	   bool verify, 
	   bool level_sync, 
	   MessageGenerator msg_gen,
	   size_t vertices_per_lock) {

#ifdef CLONE
  amplusplus::transport trans = trans_passed.clone(); // Clone transport for this run
#endif

  if (trans.rank() == 0) std::cout << "Initializing SV-CC algorithm ...  " << std::endl;

  // Instantiate algorithms used later here so we don't time their constructions
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;

  // Parent, component, and lock maps for SV CC
  std::vector<Vertex> parentS(num_vertices(g), graph_traits<Graph>::null_vertex());
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap> ParentMap;
  ParentMap parent(parentS.begin(), get(vertex_index, g));

  std::vector<int> sv_componentS(num_vertices(g), std::numeric_limits<int>::max());
  typedef iterator_property_map<std::vector<int>::iterator, VertexIndexMap> ComponentMap;
  ComponentMap component(sv_componentS.begin(), get(vertex_index, g));

  typedef boost::parallel::lock_map<VertexIndexMap> LockMap;
  LockMap locks(get(vertex_index, g), num_vertices(g) / vertices_per_lock);

  if (trans.rank() == 0) std::cout << "Number of vertices : " << num_vertices(g) << std::endl;

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(parent, v, v);
  }

  trans.set_nthreads(num_threads);
  boost::graph::distributed::sv_cc<Graph, ParentMap, MessageGenerator>
    CC(trans, g, parent, locks, msg_gen);
  trans.set_nthreads(1);

  if (level_sync) 
    CC.set_level_sync();
	  
  { amplusplus::scoped_epoch epoch(trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif


  if (trans.rank() == 0) std::cout << "Starting SV-CC algorithm ...  " << std::endl;  

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

  if (trans.rank() == 0) std::cout << "Done with SV-CC algorithm in : " << (end-start) << std::endl;

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef NUMBER_COMPONENTS
  int num_components = CC.template number_components<ComponentMap>(component);
  if (trans.rank() == 0)
    std::cout << num_components << " components found\n";
#endif

  CC.print_stats();

  if (verify) {
    if (trans.rank() == 0) std::cout << "Verifying SV-CC algorithm ...  " << std::endl;

    parent.set_consistency_model(boost::parallel::cm_forward);
    parent.set_max_ghost_cells(0);

    {
      amplusplus::scoped_epoch epoch(g.transport());
	      
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
	
	if (out_degree(v, g) == 0) {
	  assert(get(parent, v) == v);
	}

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



template <typename Graph, typename MessageGenerator>
time_type
run_sv_cc_optimized(amplusplus::transport& trans, Graph& g, 
	   int num_threads, 
	   bool verify, 
	   bool level_sync, 
	   MessageGenerator msg_gen,
	   size_t vertices_per_lock) {

#ifdef CLONE
  amplusplus::transport trans = trans_passed.clone(); // Clone transport for this run
#endif

  if (trans.rank() == 0) std::cout << "Initializing SV-CC algorithm ...  " << std::endl;

  // Instantiate algorithms used later here so we don't time their constructions
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;

  // Parent, component, and lock maps for SV CC
  std::vector<Vertex> parentS(num_vertices(g), graph_traits<Graph>::null_vertex());
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap> ParentMap;
  ParentMap parent(parentS.begin(), get(vertex_index, g));

  std::vector<int> sv_componentS(num_vertices(g), std::numeric_limits<int>::max());
  typedef iterator_property_map<std::vector<int>::iterator, VertexIndexMap> ComponentMap;
  ComponentMap component(sv_componentS.begin(), get(vertex_index, g));

  typedef boost::parallel::lock_map<VertexIndexMap> LockMap;
  LockMap locks(get(vertex_index, g), num_vertices(g) / vertices_per_lock);

  if (trans.rank() == 0) std::cout << "Number of vertices : " << num_vertices(g) << std::endl;

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(parent, v, v);
  }

  trans.set_nthreads(num_threads);
  boost::graph::distributed::connected_components<Graph, ParentMap, MessageGenerator>
    CC(trans, g, parent, locks, msg_gen);
  trans.set_nthreads(1);

  if (level_sync) 
    CC.set_level_sync();
	  
  { amplusplus::scoped_epoch epoch(trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
  epoch_times.clear();
  clear_buffer_stats();
#endif


  if (trans.rank() == 0) std::cout << "Starting SV-CC-Opt algorithm ...  " << std::endl;  

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

  if (trans.rank() == 0) std::cout << "Done with SV-CC algorithm in : " << (end-start) << std::endl;

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif

#ifdef NUMBER_COMPONENTS
  int num_components = CC.template number_components<ComponentMap>(component);
  if (trans.rank() == 0)
    std::cout << num_components << " components found\n";
#endif
  
  if (verify) {
    if (trans.rank() == 0) std::cout << "Verifying SV-CC algorithm ...  " << std::endl;

    parent.set_consistency_model(boost::parallel::cm_forward);
    parent.set_max_ghost_cells(0);

    {
      amplusplus::scoped_epoch epoch(g.transport());
	      
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
	
	if (out_degree(v, g) == 0) {
	  assert(get(parent, v) == v);
	}

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

template <typename PriorityQueueGenerator = boost::graph::distributed::cc_default_priority_queue_gen, 
	  typename Graph, 
	  typename MessageGenerator>
time_type
run_data_driven_cc(amplusplus::transport& trans, 
		   amplusplus::transport& barrier_trans, 
		   Graph& g, 
		   int num_threads, 
		   typename graph_traits<Graph>::vertices_size_type n, 
		   bool verify, 
		   bool print_stats,
		   bool allstats,
		   MessageGenerator msg_gen, 
		   MessageGenerator priority_msg_gen, 
		   unsigned int flushFreq, 
		   unsigned int eager_limit) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  if (trans.rank() == 0) std::cout << "Initializing data driven connected components ..." << std::endl;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

  std::vector<Vertex> componentS(num_vertices(g), std::numeric_limits<Vertex>::max());
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap>  ComponentsMap;
  ComponentsMap components(componentS.begin(), get(vertex_index, g));

  //#ifndef CC_PRIORITY
  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(components, v, v);
  }
  //#endif

  std::cout << "Rank: " << trans.rank() << ", Number of vertices : " << num_vertices(g)
	    << "Number of edges : " << num_edges(g) << " total vertices : " << n << std::endl;  

  trans.set_nthreads(num_threads);
  boost::graph::distributed::data_driven_cc<Graph, ComponentsMap, PriorityQueueGenerator, MessageGenerator>
    D(g, components, trans, msg_gen, priority_msg_gen, flushFreq, eager_limit);
  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  
  if (trans.rank() == 0) std::cout << "Invoking DD-CC algorithm ......" << std::endl;

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

  time_type start = D.get_start_time();  
  time_type end = D.get_end_time();

  if (trans.rank() == 0) std::cout << "Done with DD-CC algorithm in : " << (end-start) << std::endl;  

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
#endif
  
#ifdef PBGL2_PRINT_WORK_STATS
  //work_stats_t work_stats = D.get_work_stats();
  //print_and_accumulate_work_stats(work_stats, dc_stats);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    verify_cc(g, components);
  }

  if (print_stats)
    calculate_cc_stats(g, components, n, allstats);

  if (trans.rank() == 0) std::cout << "End invoking algorithm ..." << std::endl;
  return end - start;
}


// Run chaotic CC
template <typename Graph, typename MessageGenerator>
time_type
run_chaotic_cc(amplusplus::transport& trans, 
	       amplusplus::transport& barrier_trans, 
	       Graph& g, 
	       int num_threads, 
	       typename graph_traits<Graph>::vertices_size_type n, 
	       bool verify, 
	       bool print_stats,
	       bool allstats,
	       MessageGenerator msg_gen) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  if (trans.rank() == 0) std::cout << "Initializing data driven connected components ..." << std::endl;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

  std::vector<Vertex> componentS(num_vertices(g), std::numeric_limits<Vertex>::max());
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap>  ComponentsMap;
  ComponentsMap components(componentS.begin(), get(vertex_index, g));


  if (trans.rank() == 0) std::cout << "Number of vertices : " << num_vertices(g) << " total : " << n << std::endl;  

  trans.set_nthreads(num_threads);
  boost::graph::distributed::chaotic_cc<Graph, ComponentsMap, MessageGenerator>
    D(g, components, trans, n, msg_gen);
  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  
  if (trans.rank() == 0) std::cout << "Invoking DD-CC algorithm ......" << std::endl;
  time_type start = get_time();

  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
  
  D(0);
  
  // TODO: Add exception handling.  Excptions should be caught in every thread and checked in the main thread.
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();

  if (trans.rank() == 0) std::cout << "Done with Chaotic-CC algorithm in : " << (end-start) << std::endl;  

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    verify_cc(g, components);
  } // end verify

  if (print_stats)
    calculate_cc_stats(g, components, n, allstats);


  if (trans.rank() == 0) std::cout << "End invoking algorithm ..." << std::endl;
  return end - start;
}

/**
* One-to-one mapping of vertices. In this case
* vertex is mapped to running id. Id is increase column wise.
* 1  4
* 2  5
* 3  6
*/
template<typename Graph>
struct block_id_distribution {
  typedef typename graph_traits<Graph>::vertices_size_type VerticesSzType;
public:
  block_id_distribution(Graph& _pg, VerticesSzType& pn):g(_pg), n(pn) {
    if (_RANK == 0)
      std::cout << "[INFO] Using vertical id distribution" << std::endl;
  }

  
  template<typename SizeType>
  SizeType operator()(SizeType k) {
    SizeType blocksz = g.distribution().block_size(0, n);
    return (get(get(vertex_owner, g), k))*blocksz + g.distribution().local(k);
  }

private:
  Graph& g;
  VerticesSzType n;
};


/**
* One-to-one mapping of vertices. In this case
* vertex is mapped to running id. Id is increase row wise.
* 1  2
* 3  4
* 5  6
*/
template<typename Graph>
struct row_id_distribution {
public:
  row_id_distribution(Graph& _pg, int ranks):g(_pg), totalranks(ranks) {
    if (_RANK == 0)
      std::cout << "[INFO] Using horizontal id distribution" << std::endl;
  }

  template<typename SizeType>
  SizeType operator()(SizeType k) {
    auto offset = g.distribution().local(k) * totalranks;
    return offset + get(get(vertex_owner, g), k);
  }

private:
  Graph& g;
  int totalranks;
};

// Run delta-stepping CC
template <typename Graph, typename MessageGenerator, typename IdDistribution>
time_type
run_delta_cc(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g, 
	     int num_threads, 
	     typename graph_traits<Graph>::vertices_size_type n,
	     size_t nbuckets, 
	     IdDistribution idd,
	     bool verify, 
	     bool print_stats,
	     bool allstats,
	     MessageGenerator msg_gen) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  if (trans.rank() == 0) 
    std::cout << "Initializing delta connected components with buckets = " 
				   << nbuckets << "..." << std::endl;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

#ifdef PRINT_DEBUG

  block_id_distribution<Graph> logical_block_id(g, n);
  row_id_distribution<Graph> logical_row_block_id(g, trans.size());

  if (trans.rank()==0) {
    std::cout << "Printing vertex information in rank " << trans.rank() << std::endl;
    std::cout << "==========================================================" << std::endl;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      std::cout << "vid: " << v << ", local id : " << g.distribution().local(v) 
		<< ", logical block id: " << logical_block_id(v)
		<< ", logical row block id : " << logical_row_block_id(v)
		<< std::endl;

      /*      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	Vertex u = target(e, g);
	std::cout << "neighbor vid: " << u << ", local id : " << g.distribution().local(u) 
		  << ", logical block id: " << logical_block_id(u, g, n)
		  << ", logical row block id : " << logical_row_block_id(u, g, trans.size())
		  << std::endl;

		  }*/

    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (trans.rank()==1) {
    std::cout << "Printing vertex information in rank " << trans.rank() << std::endl;
    std::cout << "==========================================================" << std::endl;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      std::cout << "vid: " << v << ", local id : " << g.distribution().local(v) 
		<< ", logical block id: " << logical_block_id(v)
		<< ", logical row block id : " << logical_row_block_id(v)
		<< std::endl;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (trans.rank()==2) {
    std::cout << "Printing vertex information in rank " << trans.rank() << std::endl;
    std::cout << "==========================================================" << std::endl;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      std::cout << "vid: " << v << ", local id : " << g.distribution().local(v) 
		<< ", logical block id: " << logical_block_id(v, g, n)
		<< ", logical row block id : " << logical_row_block_id(v, g, trans.size())
		<< std::endl;

      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	Vertex u = target(e, g);
	std::cout << "neighbor vid: " << u << ", local id : " << g.distribution().local(u) 
		  << ", logical block id: " << logical_block_id(u, g, n)
		  << ", logical row block id : " << logical_row_block_id(u, g, trans.size())
		  << std::endl;

      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (trans.rank()==3) {
    std::cout << "Printing vertex information in rank " << trans.rank() << std::endl;
    std::cout << "==========================================================" << std::endl;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      std::cout << "vid: " << v << ", local id" << g.distribution().local(v) 
		<< ", logical block id: " << logical_block_id(v, g, n)
		<< ", logical row block id : " << logical_row_block_id(v, g, trans.size())
		<< std::endl;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif


  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

  std::vector<Vertex> componentS(num_vertices(g), std::numeric_limits<Vertex>::max());
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap>  ComponentsMap;
  ComponentsMap components(componentS.begin(), get(vertex_index, g));

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(components, v, v);
  }

  if (trans.rank() == 0) std::cout << "Number of vertices : " << num_vertices(g) << " total : " << n << std::endl;  

  trans.set_nthreads(num_threads);
  boost::graph::distributed::delta_stepping_cc<Graph, 
					       ComponentsMap, 
					       IdDistribution,
					       MessageGenerator>
    D(g, components, trans, nbuckets, idd,  msg_gen);

  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  
  if (trans.rank() == 0) std::cout << "Invoking Delta-CC algorithm ......" << std::endl;
  time_type start = get_time();

  // Many threads now
  trans.set_nthreads(num_threads);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (int i = 0; i < num_threads - 1; ++i) {
    boost::thread thr(boost::ref(D), i + 1);
    threads[i].swap(thr);
  }
  
  D(0);
  
  // TODO: Add exception handling.  Excptions should be caught in every thread and checked in the main thread.
  
  for (int i = 0; i < num_threads - 1; ++i)
    threads[i].join();
  
  time_type end = get_time();

  if (trans.rank() == 0) std::cout << "Done with Delta-CC algorithm in : " << (end-start) << std::endl;  

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    verify_cc(g, components);
  }

  if (print_stats)
    calculate_cc_stats(g, components, n, allstats);

  if (trans.rank() == 0) std::cout << "End invoking algorithm ..." << std::endl;
  return end - start;
}


// Run level CC
template <typename Graph, typename MessageGenerator>
time_type
run_level_sync_cc(amplusplus::transport& trans, 
		  amplusplus::transport& barrier_trans, 
		  Graph& g, 
		  int num_threads, 
		  typename graph_traits<Graph>::vertices_size_type n,
		  bool verify, 
		  bool print_stats,
		  bool allstats,
		  MessageGenerator msg_gen) {
#ifdef CLONE
  amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif
  if (trans.rank() == 0) 
    std::cout << "Initializing level sync connected components... " 
				   << std::endl;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

  std::vector<Vertex> componentS(num_vertices(g), std::numeric_limits<Vertex>::max());
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap>  ComponentsMap;
  ComponentsMap components(componentS.begin(), get(vertex_index, g));

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    put(components, v, v);
  }

  if (trans.rank() == 0) std::cout << "Number of vertices : " << num_vertices(g) << " total : " << n << std::endl;  

  trans.set_nthreads(num_threads);
  boost::graph::distributed::level_sync_cc<Graph,
					   ComponentsMap, 
					   MessageGenerator>
    D(g, components, trans, n,  msg_gen);

  trans.set_nthreads(1);

  { amplusplus::scoped_epoch epoch(barrier_trans); }
  
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  epoch_times.clear();
  clear_buffer_stats();
#endif
  
  if (trans.rank() == 0) std::cout << "Invoking level-sync-CC algorithm ......" << std::endl;

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
  time_type start = D.get_start_time();

  if (trans.rank() == 0) 
    std::cout << "Done with level-sync-CC algorithm in : " << (end-start) 
	      << ", buckets(levels) processed : " 
	      << D.get_num_levels() << std::endl;  

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
  print_and_clear_epoch_times();
  print_buffer_stats();
#endif
  
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  //const std::pair<unsigned long long, unsigned long long> stats = D.get_cache_stats();
  //print_and_accumulate_cache_stats(stats);
#endif
  
  // Back to one thread
  trans.set_nthreads(1);

  if (verify) {
    verify_cc(g, components);
  }

  if (print_stats)
    calculate_cc_stats(g, components, n, allstats);

  if (trans.rank() == 0) std::cout << "End invoking algorithm ..." << std::endl;
  return end - start;
}



enum mode_type {mode_none, mode_self, mode_async_bfs, mode_ds_async_bfs, mode_level_synchronized_bfs, mode_delta_stepping, mode_dc, mode_connected_components, mode_page_rank};
enum routing_type {rt_none, rt_hypercube, rt_rook};

// functions to separate edge weight and edge

// separating edge
template<typename t1, typename t2>
struct select_edge : std::unary_function< std::pair<t1, std::pair<t1, t2> >, std::pair<t1, t1> > {
public:
  select_edge() {}

  std::pair<t1, t1> operator()(std::pair<t1, std::pair<t1, t2> >& ew) const {
    return std::make_pair(ew.first, ew.second.first);
  }
};

// separating edge weight
template<typename t1, typename t2>
struct select_edge_weight : std::unary_function< std::pair<t1, std::pair<t1, t2> >, t2 > {
public:
  select_edge_weight() {}

  t2 operator()(std::pair<t1, std::pair<t1, t2> >& ew) const {
    return ew.second.second;
  }
};

template<typename T>
void time_statistics(std::vector<T>& data, 
		     T& m, 
		     T& minimum,
		     T& q1,
		     T& median,
		     T& q3,
		     T& maximum,
		     T& stddev) {

  using namespace boost::accumulators;
  accumulator_set<T, stats<tag::variance> > acc;
  for_each(data.begin(), data.end(), bind<void>(ref(acc), _1));
  
  m = mean(acc);
  stddev = sqrt(variance(acc));

  auto const min = 0;
  auto const Q1 = data.size() / 4;
  auto const Q2 = data.size() / 2;
  auto const Q3 = Q1 + Q2;
  auto const max = data.size() - 1;


  std::nth_element(data.begin(),          data.begin() + min, data.end());
  minimum = data[min];
  std::nth_element(data.begin(),          data.begin() + Q1, data.end());
  q1 = data[Q1];
  std::nth_element(data.begin() + Q1 + 1, data.begin() + Q2, data.end());
  median = data[Q2];
  std::nth_element(data.begin() + Q2 + 1, data.begin() + Q3, data.end()); 
  q3 = data[Q3];
  std::nth_element(data.begin() + Q3 + 1, data.begin() + max, data.end()); 
  maximum = data[max];

}
		     

class dc_test {

private:
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
  std::vector<routing_type> routing;
  double edge_list_reserve_factor = 1.15;
  bool no_reductions = false, per_thread_reductions = true;
  std::vector<size_t> coalescing_size;
  size_t distribution_coalescing_size = 1 << 17;
  std::vector<size_t> reduction_cache_size;
  std::vector<unsigned int> number_poll_task;
  std::vector<std::string> dc_data_structures;
  std::vector<unsigned int> flushFreq;
  std::vector<unsigned int> eager_limit;
  std::vector<unsigned int> recvDepth;
  std::vector<unsigned int> delta;
  std::vector<unsigned int> k_levels = {1};
  std::vector<unsigned int> priority_coalescing_size_v = {40000};
  bool run_dc = false, run_chaotic = false, run_ds = false, run_kla = false;

  // distributions
  enum distribution_t {
    block,
    cyclic
  };

  size_t block_size = 10;

  // default block distribution
  distribution_t distribution_type = block;
  id_distribution_t id_distribution;

  // mis
  bool run_ss_mis = false;
  bool run_bucket_mis = false;
  bool run_luby_mis = false;
  std::vector<std::string> luby_algorithms;
  
  // connected components
  bool run_cc = false;
  std::vector<std::string> cc_algorithms;
  // Available CC algorithms
  // 1. run_sv_cc
  // 2. run_sv_cc_level_sync
  // 3. run_sv_cc_opt
  // 4. run_sv_cc_opt_level_sync
  // 5. run_sv_ps_cc // optimized version
  // 6. run_sv_ps_cc_level_sync // optmized version
  // 7. run_dd_cc
  // 8. run_cc_chaotic
  // 9. run_cc_ds // delta stepping cc
  // 10. run_level_sync_cc // level sync agm
  // data driven connected components
  size_t vertices_per_lock = 64;
  std::vector<size_t> nbuckets;
  size_t cutoff_degree = 0;
  // CC stats
  // print CC stats
  bool print_cc_stats = false;
  // print all CC stats (otherwise limited to 5)
  bool allstats = false;

  // read graph from file
  std::string graph_file;
  bool read_graph = false;
  bool gen_graph = false;
  bool with_weight = false;


public:
  dc_test(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if(arg == "--receive-depth") {
	recvDepth = extract_params<unsigned int>( argv[i+1] );
	//std::cerr<<"receiveDepth from main"<<recvDepth<<std::endl;
      }
      if (arg == "--threads") {
	thread_num_vals = extract_params<int>(argv[i+1]);
      }

      if (arg == "--graph-file") {
	graph_file = argv[i+1];
	read_graph = true;
      }

      if (arg == "--with-weight") {
	// read the graph file with weight
	with_weight = true;
      }

      if (arg == "--seed") {
	seed64 = boost::lexical_cast<uint64_t>( argv[i+1] );
      }

      if (arg == "--scale") {
	scale = boost::lexical_cast<size_t>( argv[i+1] );
	n = (unsigned long long)(1) << scale;
	gen_graph = true;	
      }

      if (arg == "--degree")
	edgefactor = boost::lexical_cast<double>( argv[i+1] );

      if (arg == "--num-sources")
	num_sources = boost::lexical_cast<size_t>( argv[i+1] );

      if (arg == "--iterations")
	iterations = boost::lexical_cast<size_t>( argv[i+1] );

      if (arg == "--verify")
	verify = true;

      if (arg == "--stats")
	stats = true;

      if (arg == "--coalescing-size") {
	coalescing_size = extract_params<size_t>( argv[i+1] );
      }

      if (arg == "--reduction-cache-size")
        reduction_cache_size = extract_params<size_t>( argv[i+1] );

      if (arg == "--distribution-coalescing-size")
	distribution_coalescing_size = boost::lexical_cast<size_t>( argv[i+1] );

      if (arg == "--distribution") {
	if (strcmp(argv[i+1],"block") == 0)
	  distribution_type = block;
	else if (strcmp(argv[i+1],"cyclic") == 0)
	  distribution_type = cyclic;
	else {
	  std::cout << "Invalid distribution type. Available types are block and cyclic" << std::endl;
	  //	  abort = true;	  
	}
      }

      if (arg == "--block-size") {
	block_size = boost::lexical_cast<size_t>( argv[i+1] );
      }

      if (arg == "--id-distribution") {
	if (strcmp(argv[i+1],"vertical") == 0)
	  id_distribution = vertical;
	else if (strcmp(argv[i+1],"horizontal") == 0)
	  id_distribution = horizontal;
	else {
	  std::cout << "Invalid id distribution type. Available types are vertical and horizontal" 
		    << std::endl;
	  //	  abort = true;	  
	}
      }

      if (arg == "--max-weight")
	C = boost::lexical_cast<weight_type>( argv[i+1] );

      if(arg == "--poll-task"){
	number_poll_task = extract_params<unsigned int>( argv[i+1] );
	//std::cerr<<"poll task number:" <<number_poll_task<<"\n";
      }

      if (arg == "--flush") {
	flushFreq = extract_params<unsigned int> ( argv[i+1] );
      }

      if (arg == "--eager-limit") {
	eager_limit = extract_params<unsigned int> ( argv[i+1] );
      }

      if (arg == "--priority_coalescing_size") {
	priority_coalescing_size_v = extract_params<unsigned int> ( argv[i+1] );
      }

      if(arg == "--ds"){
	dc_data_structures = extract_params<std::string> ( argv[i+1] );
      }

      if (arg == "--rook")
	routing.push_back(rt_rook);

      if (arg == "--rt_none")
	routing.push_back(rt_none);

      if (arg == "--with-no-reductions")
	no_reductions = true;

      if (arg == "--without-no-reductions")
	no_reductions = false;

      if (arg == "--with-per-thread-reductions")
	per_thread_reductions = true;

      if (arg == "--without-per-thread-reductions")
	per_thread_reductions = false;

      if (arg == "--run_dc")
	run_dc = true;

      if (arg == "--run_chaotic")
	run_chaotic = true;

      // CC
      if (arg == "--run_cc")
	run_cc = true;

      if (arg == "--cc_algorithms") {
	cc_algorithms = extract_params<std::string> ( argv[i+1] );
      }

      if (arg == "--vertices-per-lock") {
	if (!run_cc) {
	  std::cerr << "Vertices per lock is only allowed for Shiloach-Vishky connected components." 
		    << std::endl;
	}
	vertices_per_lock = boost::lexical_cast<size_t>( argv[i+1] );
      }

      if (arg == "--level-sync") {
	if (!run_cc) {
	  std::cerr << "Level synch is only allowed for Shiloach-Vishky connected components." 
		    << std::endl;
	  return;
	}
	level_sync = true;
      }

      if (arg == "--buckets") {
	nbuckets = extract_params<size_t> ( argv[i+1] );
      }

      if (arg == "--print_cc_stats")
	print_cc_stats = true;

      if (arg == "--allstats")
	allstats = true;

      if (arg == "--cutoff_degree") {
	cutoff_degree = boost::lexical_cast<size_t>( argv[i+1] );
      }

      // MIS
      if (arg == "--run_ss_mis")
	run_ss_mis = true;

      if (arg == "--run_bucket_mis")
	run_bucket_mis = true;

      if (arg == "--run_luby_mis") {
	run_luby_mis = true;
      }

      if (arg == "--luby_algorithms") {
	luby_algorithms = extract_params<std::string> ( argv[i+1] );
      }

      if (arg == "--run_ds")
	run_ds = true;
      
      if (arg == "--run_kla")
	run_kla = true;
      
      if (arg == "--klevel")
	k_levels = extract_params<unsigned int> ( argv[i+1] );

      if (arg == "--delta") {
	delta = extract_params<unsigned int> ( argv[i+1] );
      }
    }


    if (thread_num_vals.empty()) thread_num_vals.push_back(1);
    if (routing.empty()) routing.push_back(rt_none);

    // Check for necessary parameters
    bool abort = false;

    if (read_graph && gen_graph) {
      std::cout << "Both generate graph and read graph from file options are specified. Aborting !"
		<< std::endl;
      abort = true;
    }

    if(coalescing_size.empty()) {
      std::cerr << "Provide a list of coalescing sizes: --coalescing-size x,y,..." << std::endl;
      abort = true;
    }
    if(recvDepth.empty()) {
      std::cerr << "Provide a list of receive depths: --receive-depth x,y,..." << std::endl;
      abort = true;
    }
    if(number_poll_task.empty()) {
      std::cerr << "Provide a list of poll tasks: --poll-task x,y,..." << std::endl;
      abort = true;
    }
    if(flushFreq.empty()) {
      std::cerr << "Provide a list of flush frequencies: --flush x,y,..." << std::endl;
      abort = true;
    }
    if(!(run_dc || run_ds || run_kla || run_cc || run_chaotic || 
	 run_ss_mis || run_luby_mis || run_bucket_mis)) {
      std::cerr << "Select at least one algorithm (--run_dc or --run_ds or" 
		<< " --run_kla or --run_cc or --run_chaotic or"
                << " --run_ss_mis or --run_luby_mis or --run_bucket_mis)." << std::endl;
      abort = true;
    } else if (run_cc) {
      if (cc_algorithms.empty()) {
	std::cerr << "Must specify atleast one connected component algorithms to run. Options : "
		  << "run_sv_cc, run_sv_ps_cc, run_sv_cc_level_sync, run_sv_ps_cc_level_sync, run_dd_cc, run_cc_chaotic, run_cc_ds, run_level_sync_cc"
		  << ", run_sv_cc_opt, run_sv_cc_opt_level_sync"
		  << std::endl;
	abort = true;
      }
    }

    if (run_luby_mis && luby_algorithms.empty()) {
      std::cerr << "Select at least one Luby algorithm using --luby_algorithms. "
		<< "Available algorithms are A, AV1,AV2, B"
		<< std::endl;
      abort = true;
    }
 
    if(abort) std::abort();
  }

  typedef unsigned long long block_node_t;

  void test_main(int argc, char* argv[]) {

    amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
    amplusplus::transport trans = env.create_transport();
    amplusplus::transport barrier_trans = trans.clone();

    _RANK = trans.rank(); // For performance counter output

    if (_RANK == 0) {
#ifdef NEW_GRAPH500_SPEC
      std::cout << "Synthetic graph generator used -- RMAT-2" << std::endl;
#else
      std::cout << "Synthetic graph generator used -- RMAT-1" << std::endl;
#endif
    }


    typedef amplusplus::transport::rank_type rank_type;
    typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge, no_property, distributedS<unsigned long long> > Digraph;

    typedef graph_traits<Digraph>::vertex_descriptor Vertex;
    typedef graph_traits<Digraph>::edges_size_type edges_size_type;
    typedef graph_traits<Digraph>::vertices_size_type vertices_size_type;

    typedef property_map<Digraph, vertex_index_t>::type VertexIndexMap;
    typedef property_map<Digraph, vertex_owner_t>::const_type OwnerMap;

    // Output of permutation
    // edge is expressed as source_vertex, target_vertex
    // But in the distribution we need source vertex as first element in the pair
    // Therefore we are putting source vertex first and then another pair with target and weight
    // So it will be pair<source, pair<target, weight> >.
    typedef std::vector<std::pair<edges_size_type, std::pair<edges_size_type, weight_type> > > edge_with_weight_t;
    edge_with_weight_t edges;
    parallel::variant_distribution<vertices_size_type> distrib;

    time_type gen_start;
    // Seed general-purpose RNG
    rand48 gen, synch_gen;
    gen.seed(seed64);
    synch_gen.seed(seed64);

    typedef generator_iterator<uniform_int<weight_type>, rand48> weight_iterator_t;
    weight_iterator_t gi 
      = make_generator_iterator(gen, uniform_int<weight_type>(1, C));

    if (read_graph) {
      gen_start = get_time();
      std::cout << "Reading graph from file : " << graph_file << std::endl;
      graph_reader<block_node_t> gr(graph_file.c_str());
      gr.read_header();

      n = gr.get_num_vertices();
      edges_size_type m = gr.get_num_edges();

      // reserve space for edges
      std::cout << "Reserving space for edges : " << m << std::endl;
      edges.reserve(m);
      std::cout << "Reserved space for edges." << std::endl;

      auto bdist = parallel::block<vertices_size_type>(trans, n);
      // starting position for current rank
      block_node_t start_index = bdist.start(trans.rank());
      // locally stored ids in current rank
      block_node_t local_count = bdist.block_size(n);
      distrib = bdist;

      if (with_weight) {
	if (_RANK == 0)
	  std::cout << "Reading edges with weights ..." << std::endl;

	if (!gr.read_edges_with_weight<edges_size_type, weight_type,
	    edge_with_weight_t>(start_index, local_count, edges))
	  return;
      } else {
	if (_RANK == 0)
	  std::cout << "Reading edges with generated weights ..." << std::endl;

	if (!gr.read_edges_wo_weight<edges_size_type,edge_with_weight_t,
	    weight_iterator_t>(start_index,
			       local_count,
			       edges,
			       gi))
	  return;
      }

      std::cout << "Total number of edges read : " << edges.size()
		<< std::endl;


    } else {
      assert(gen_graph);

      edges_size_type m = static_cast<edges_size_type>(floor(n * edgefactor));

      gen_start = get_time();

      edges.reserve(static_cast<edges_size_type>(floor(edge_list_reserve_factor * 2 * m / trans.size())));

      if (distribution_type == cyclic) {
	distrib = parallel::oned_block_cyclic<vertices_size_type>(trans, block_size);
	std::cout << "Distribution is cyclic. block size : " << block_size << std::endl;
      } else //(distribution_type == block)
	distrib = parallel::block<vertices_size_type>(trans, n);
    

      {
	boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());

#ifdef CLONE
	amplusplus::transport trans = trans.clone(); // Clone transport for distribution
#endif

	edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
	edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

	// Permute and redistribute copy constructs the input iterator
	uint64_t a = rand_64(gen);
	uint64_t b = rand_64(gen);

	// Build a graph to test with
	typedef graph500_iterator<Digraph, generator_iterator<uniform_int<weight_type>, rand48>, weight_type > Graph500Iter;

	// Select the highest thread val
	int num_threads = 1;
	for(unsigned int thread_choice : thread_num_vals) {
	  if (num_threads < thread_choice)
	    num_threads = thread_choice;
	}


	// As for now, we assume that only the first routing from the list is used to generate the graph.  It seems that it does not really matter which routing we use here since generation of the graph is not a part of the performance test.
	if (routing[0] == rt_none) {
	  do_distribute<Graph500Iter, amplusplus::no_routing>(distrib, trans, e_start, e_count, a, b, edges, gi, num_threads);
	} else if (routing[0] == rt_hypercube) {
	  do_distribute<Graph500Iter, amplusplus::hypercube_routing>(distrib, trans, e_start, e_count, a, b, edges, gi, num_threads);
	} else if (routing[0] == rt_rook) {
	  do_distribute<Graph500Iter, amplusplus::rook_routing>(distrib, trans, e_start, e_count, a, b, edges, gi, num_threads);
	}
      }
    }


    typedef select_edge<edges_size_type, weight_type> EdgeSelectFunction;
    typedef transform_iterator<EdgeSelectFunction, edge_with_weight_t::iterator> edge_only_iterator;

    typedef select_edge_weight<edges_size_type, weight_type> WeightSelectFunction;
    typedef transform_iterator<WeightSelectFunction, edge_with_weight_t::iterator> weight_only_iterator;

    edge_only_iterator edge_begin(edges.begin(), EdgeSelectFunction()), 
      edge_end(edges.end(), EdgeSelectFunction());

    weight_only_iterator weight_begin(edges.begin(), WeightSelectFunction());

    std::cout << "rank -- " << trans.rank() << ", total edges received=" << edges.size() << std::endl; 
    time_type gt1 = get_time();
    Digraph g(edges_are_unsorted_multi_pass, edge_begin, edge_end,
	      weight_begin, n, trans, distrib);

    //Digraph g(edges_are_sorted, edge_begin, edge_end,
    //	      weight_begin, n, trans, distrib);
    time_type gt2 = get_time();
    std::cout << "Local graph creation time : " << (gt2-gt1)
	      << std::endl;
    // Clear edge array above
    edges.clear();

    //============= run some stats ===================//
    run_grah_stats(g);

    //Generate sources
    boost::uniform_int<uint64_t> rand_vertex(0, n-1);

    { amplusplus::scoped_epoch epoch(trans); }

    time_type gen_end = get_time();

    if (trans.rank() == 0) {
      std::cout << "Graph generation took " << print_time(gen_end - gen_start) << "s\n";

#ifdef PRINT_DEBUG
      // Printing graph edges and edge weights
      // Property maps
      typedef property_map<Digraph, weight_type WeightedEdge::*>::type WeightMap;
      WeightMap weight = get(&WeightedEdge::weight, g);

      typename graph_traits < Digraph >::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
	weight_type we = get(weight, *ei);
	Vertex v = source(*ei, g);
	std::cout << " Weight - " << we << " edge source : " << boost::source(*ei, g)
		  << " edge target : " << boost::target(*ei, g) << std::endl;
      }
#endif
      
    }

    // Max degree computation and printout
    /*{
      typedef typename graph_traits<Digraph>::degree_size_type Degree;

      // Compute the maximum edge degree
      Degree max_degree = 0;
      BGL_FORALL_VERTICES_T(u, g, Digraph) {
      max_degree = max BOOST_PREVENT_MACRO_SUBSTITUTION (max_degree, out_degree(u, g));
      }
      // max_degree = all_reduce(process_group(g), max_degree, maximum<Degree>());
      all_reduce<Degree, maximum<Degree> > r(trans, maximum<Degree>());
      max_degree = r(max_degree);

      if(trans.rank() == 0)
      std::cout << "Maximum Degree is " << max_degree << std::endl;
      }*/
      
    // Property maps
    typedef property_map<Digraph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = get(&WeightedEdge::weight, g);

    // Distance map
    std::vector<weight_type> distanceS(num_vertices(g), std::numeric_limits<weight_type>::max());
    typedef iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap>  DistanceMap;
    DistanceMap distance(distanceS.begin(), get(vertex_index, g));

    std::string empty_ds = "";

    for(unsigned int thread_choice : thread_num_vals) {
      for(unsigned int coalescing_choice : coalescing_size) {
	for(unsigned int depth_choice : recvDepth) {
	  for(unsigned int poll_choice : number_poll_task) {
	    for(routing_type routing_choice : routing) {
	      // First print the configuration
	      if(run_dc) {
		for(unsigned int flush_choice : flushFreq) {
		  for(unsigned int eager_choice : eager_limit) {
		    for(unsigned int priority_coalescing_size : priority_coalescing_size_v) {
		      for (std::string dc_ds : dc_data_structures) {
			if (trans.rank() == 0) std::cout << "[DC]" 
							 << " Data Structure :" << dc_ds
							 << " Threads: " << thread_choice 
							 << " Coalescing: " << coalescing_choice 
							 << " Poll: " << poll_choice 
							 << " Routing: " << routing_choice 
							 << " Depth: " << depth_choice;
			select_alg(env, 
				   thread_choice, 
				   coalescing_choice, 
				   priority_coalescing_size, 
				   depth_choice, 
				   poll_choice, 
				   routing_choice, 
				   g, 
				   distance, 
				   weight, 
				   thread_choice, 
				   flush_choice, 
				   eager_choice, 
				   0, 
				   dc_ds, 
				   synch_gen, 
				   rand_vertex, 
				   false, //kla
				   true, //dc
				   false, //ds
				   false, //cc
				   false, //chaotic
				   false, //ss-mis
				   false, // bucket mis
				   false, //luby-mis
				   (size_t)1);
		      }
		    }
		  }
		}
	      }
	      if (run_ss_mis || run_luby_mis || run_bucket_mis) {
		for(unsigned int flush_choice : flushFreq) {
		  for (std::string dc_ds : dc_data_structures) {
		    if (run_ss_mis) {
		      if (trans.rank() == 0) std::cout << "[MIS-SS]" 
						       << " Data Structure :" << dc_ds
						       << " Threads: " << thread_choice 
						       << " Coalescing: " << coalescing_choice 
						       << " Poll: " << poll_choice 
						       << " Routing: " << routing_choice 
						       << " Depth: " << depth_choice
						       << " Flush: " << flush_choice;
		      
		      select_alg(env, 
				 thread_choice, 
				 coalescing_choice, 
				 10000, 
				 depth_choice, 
				 poll_choice, 
				 routing_choice, 
				 g, 
				 distance, 
				 weight, 
				 thread_choice, 
				 flush_choice, 
				 0, 
				 0, 
				 dc_ds, 
				 synch_gen, 
				 rand_vertex, 
				 false, //kla
				 false,  // dc
				 false, // ds 
				 false, // cc
				 false, // chaotic
				 true, // ss-mis
				 false, // bucket-mis
				 false, // luby-mis
				 1);
		    }

		    if (run_bucket_mis) {
		      if (trans.rank() == 0) std::cout << "[MIS-BUCKET]" 
						       << " Data Structure :" << dc_ds
						       << " Threads: " << thread_choice 
						       << " Coalescing: " << coalescing_choice 
						       << " Poll: " << poll_choice 
						       << " Routing: " << routing_choice 
						       << " Depth: " << depth_choice
						       << " Flush: " << flush_choice;
		      
		      select_alg(env, 
				 thread_choice, 
				 coalescing_choice, 
				 10000, 
				 depth_choice, 
				 poll_choice, 
				 routing_choice, 
				 g, 
				 distance, 
				 weight, 
				 thread_choice, 
				 flush_choice, 
				 0, 
				 0, 
				 dc_ds, 
				 synch_gen, 
				 rand_vertex, 
				 false, //kla
				 false,  // dc
				 false, // ds 
				 false, // cc
				 false, // chaotic
				 false, // ss-mis
				 true, //bucket-mis
				 false, // luby-mis
				 1);
		    }


		    if (run_luby_mis) {
		      for (std::string luby_algo : luby_algorithms) {
			if (trans.rank() == 0) std::cout << "[MIS-LUBY]" 
		       					 << " Luby Algorithm :" << luby_algo
							 << " Data Structure :" << dc_ds
							 << " Threads: " << thread_choice 
							 << " Coalescing: " << coalescing_choice 
							 << " Poll: " << poll_choice 
							 << " Routing: " << routing_choice 
							 << " Depth: " << depth_choice
							 << " Flush: " << flush_choice;

			select_alg(env, 
				   thread_choice, 
				   coalescing_choice, 
				   10000, 
				   depth_choice, 
				   poll_choice, 
				   routing_choice, 
				   g, 
				   distance, 
				   weight, 
				   thread_choice, 
				   flush_choice, 
				   0, 
				   0, 
				   dc_ds, 
				   synch_gen, 
				   rand_vertex, 
				   false, //kla
				   false,  // dc
				   false, // ds 
				   false, // cc
				   false, // chaotic
				   false, // ss-mis
				   false, //bucket mis
				   true, // luby-mis
				   1,
				   luby_algo);

		      }
		    }
		  }
		}
	      }
	      if(run_chaotic) {
		for(unsigned int flush_choice : flushFreq) {
		  for(unsigned int eager_choice : eager_limit) {
		    if (trans.rank() == 0) std::cout << "[Chaotic]"
						     << " Threads: " << thread_choice 
						     << " Coalescing: " << coalescing_choice 
						     << " Poll: " << poll_choice 
						     << " Routing: " << routing_choice 
						     << " Depth: " << depth_choice;
		    select_alg(env, 
			       thread_choice, 
			       coalescing_choice, 
			       coalescing_choice, 
			       depth_choice, 
			       poll_choice, 
			       routing_choice, 
			       g, 
			       distance, 
			       weight, 
			       thread_choice, 
			       flush_choice, 
			       eager_choice, 
			       0, 
			       empty_ds, 
			       synch_gen, 
			       rand_vertex, 
			       false, //kla 
			       false,  // dc
			       false,  // ds
			       false,  // cc
			       true,  // chaotic
			       false, //ss-mis
			       false, // bucket mis
			       false, //luby-mis
			       1);
		  }
		}
	      }
	      if(run_ds) { // run_ds
		for(unsigned int delta_choice : delta) {
		  for (std::string dc_ds : dc_data_structures) {
		    for(unsigned int flush_choice : flushFreq) {
		      if (trans.rank() == 0) std::cout << "[Delta]" 
						       << " Data Structure :" << dc_ds
						       << " Threads: " << thread_choice
						       << " Flush : " << flush_choice  
						       << " Coalescing: " << coalescing_choice 
						       << " Poll: " << poll_choice 
						       << " Routing: " << routing_choice 
						       << " Depth: " << depth_choice;
		      select_alg(env, 
				 thread_choice, 
				 coalescing_choice, 
				 0, 
				 depth_choice, 
				 poll_choice, 
				 routing_choice, 
				 g, 
				 distance, 
				 weight, 
				 thread_choice, 
				 flush_choice, 
				 0, 
				 delta_choice, 
				 dc_ds, 
				 synch_gen, 
				 rand_vertex, 
				 false, //kla
				 false, //dc
				 true,  // ds
				 false, // cc
				 false, // chaotic
				 false, //ss-mis
				 false, // bucket mis
				 false, //ss-luby
				 1);
		    }
		  }
		}
	      }
	      if(run_kla) {
		unsigned int delta_choice = 1;
		for(unsigned int k_level : k_levels) {
		  for (std::string dc_ds : dc_data_structures) {
		    for(unsigned int flush_choice : flushFreq) {
		      if (trans.rank() == 0) std::cout << "[KLA]" 
						       << " Data Structure :" << dc_ds
						       << " Threads: " << thread_choice 
						       << " Flush : " << flush_choice 
						       << " Coalescing: " << coalescing_choice 
						       << " Poll: " << poll_choice 
						       << " Routing: " << routing_choice 
						       << " Depth: " << depth_choice 
						       << "k_level: " << k_level;
		      select_alg(env, 
				 thread_choice, 
				 coalescing_choice, 
				 0, 
				 depth_choice, 
				 poll_choice, 
				 routing_choice, 
				 g, 
				 distance, 
				 weight, 
				 thread_choice, 
				 flush_choice, 
				 0, 
				 delta_choice, 
				 dc_ds, 
				 synch_gen, 
				 rand_vertex, 
				 true, //kla
				 false, //dc
				 false, //ds
				 false, //cc
				 false, //chaotic
				 false, //ss-mis
				 false, // bucket mis
				 false, //luby
				 k_level);		
		    }
		  }
		}
	      }

	      if (run_cc) {
		for (std::string algorithm : cc_algorithms) {
		  if (algorithm == "run_dd_cc") {
		    for(unsigned int flush_choice : flushFreq) {
		      if (trans.rank() == 0) std::cout << "[CC]" 
						       << " Algorithm : " << algorithm 
						       << " Threads: " << thread_choice 
						       << " Coalescing: " << coalescing_choice 
						       << " Poll: " << poll_choice 
						       << " Routing: " << routing_choice 
						       << " Depth: " << depth_choice
						       << " Flush: " << flush_choice;

		      unsigned int delta_choice = 1;
		      select_alg(env, 
				 thread_choice, 
				 coalescing_choice, 
				 1, 
				 depth_choice, 
				 poll_choice, 
				 routing_choice, 
				 g, 
				 distance, 
				 weight, 
				 thread_choice, 
				 0, 
				 0, 
				 delta_choice, 
				 empty_ds, 
				 synch_gen, 
				 rand_vertex, 
				 false, //kla
				 false, //dc
				 false, //ds
				 true, //cc
				 false, //chaotic
				 false, //ss-mis
				 false, //bucket mis
				 false, //luby
				 0, // k-level
				 "", // luby algo
				 algorithm);
		    }
		  } else {

		    if (trans.rank() == 0) std::cout << "[CC]" 
						     << " Algorithm : " << algorithm 
						     << " Threads: " << thread_choice 
						     << " Coalescing: " << coalescing_choice 
						     << " Poll: " << poll_choice 
						     << " Routing: " << routing_choice 
						     << " Depth: " << depth_choice;

		    if ("run_cc_ds" == algorithm) {
		      unsigned int delta_choice = 1; // ignore
		      for(size_t bucketsz : nbuckets) {
			select_alg(env, 
				   thread_choice, 
				   coalescing_choice, 
				   1, 
				   depth_choice, 
				   poll_choice, 
				   routing_choice, 
				   g, 
				   distance, 
				   weight, 
				   thread_choice, 
				   0, 
				   0, 
				   delta_choice, 
				   empty_ds, 
				   synch_gen, 
				   rand_vertex, 
				   false, //kla
				   false, //dc
				   false, //ds
				   true, //cc
				   false, //chaotic
				   false, //ss-mis
				   false, // bucket mis
				   false, //luby
				   0, // k-level
				   "", // luby algo
				   algorithm,
				   bucketsz);	
		      }
		    } else {
		      // delta is ignored
		      unsigned int delta_choice = 1;
		      select_alg(env, 
				 thread_choice, 
				 coalescing_choice, 
				 1, 
				 depth_choice, 
				 poll_choice, 
				 routing_choice, 
				 g, 
				 distance, 
				 weight, 
				 thread_choice, 
				 0, 
				 0, 
				 delta_choice, 
				 empty_ds, 
				 synch_gen, 
				 rand_vertex, 
				 false, //kla
				 false, //dc
				 false, //ds
				 true, //cc
				 false, //chaotic
				 false, //ss-mis
				 false, // bucket mis
				 false, //luby
				 0, // k-level
				 "", // luby algo
				 algorithm);	
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

private:
  template<typename Result>
  std::vector<Result> extract_params(const std::string args) {
    size_t d = 0, d2;

    std::vector<Result> r;
    while ((d2 = args.find(',', d)) != std::string::npos) {
      r.push_back(boost::lexical_cast<Result>(args.substr(d, d2 - d)));
      d = d2 + 1;
    }

    r.push_back(boost::lexical_cast<Result>(args.substr(d, args.length())));

    return r;
  }

  template<typename Digraph, typename DistanceMap, typename WeightMap, typename RNG, typename RandVertex>
  void select_alg(amplusplus::environment& env, 
		  unsigned int thread_choice, 
		  unsigned int coalescing_choice, 
		  unsigned int priority_coalescing_size, 
		  unsigned int depth_choice, 
		  unsigned int poll_choice, 
		  routing_type routing_choice, 
		  Digraph &g, 
		  DistanceMap &distance, 
		  WeightMap &weight, 
		  unsigned int num_threads, 
		  unsigned int flush_choice, 
		  unsigned int eager_choice, 
		  unsigned int delta_choice, 
		  std::string& dc_ds_choice, 
		  RNG synch_gen, 
		  RandVertex rand_vertex, 
		  bool run_kla, 
		  bool run_dc, 
		  bool run_ds, 
		  bool run_cc, 
		  bool run_chaotic, 
		  bool run_ss_mis,
		  bool run_bucket_mis,
		  bool run_luby_mis,
		  size_t k_level,
		  std::string luby_algo="",
		  std::string cc_algo="",
		  size_t bucket_choice=0) {

    //run_ds = !run_dc;
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    flushes.reset(new atomic_flush_type[coalescing_choice]);
    flushes_size = coalescing_choice;
    all_flushes.resize(coalescing_choice);
    cumulative_flushes.resize(coalescing_choice);
#endif
    
    // We will not run two algorithms at once since they have disjoint sets of parameters
    if(run_dc && run_ds && run_kla) abort();

    amplusplus::transport trans = env.create_transport();
    amplusplus::transport barrier_trans = trans.clone();

    if(run_dc && trans.rank() == 0) {
      std::cout << " Priority: " << priority_coalescing_size << " Flush: " << flush_choice  << " Eager: " << eager_choice;
    }

    if(run_ds && trans.rank() == 0) { 
      std::cout << " Delta: " << delta_choice;
    }

    if(run_kla && trans.rank() == 0) { 
      std::cout << " k_level: " << k_level<<std::endl;
    }

    // Create transport for the algorithm
    env.template downcast_to_impl<amplusplus::detail::mpi_environment_obj>()->set_poll_tasks(poll_choice);
    env.template downcast_to_impl<amplusplus::detail::mpi_environment_obj>()->set_recv_depth(depth_choice);

    if (trans.rank() == 0)
      std::cout << "per_thread_reductions - " << per_thread_reductions << " no_reductions : " << no_reductions << std::endl;

    if(per_thread_reductions) {
      for(unsigned int reduction_choice : reduction_cache_size) {
	if (trans.rank() == 0) std::cout << " Reduction: " << reduction_choice << std::endl;
	if (routing_choice == rt_none) {
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_choice), reduction_choice, amplusplus::no_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_choice, amplusplus::no_routing(trans.rank(), trans.size()));
	  select_msg_gen(msg_gen, priority_msg_gen, trans, barrier_trans, g, distance, weight, thread_choice, flush_choice, eager_choice, delta_choice, dc_ds_choice, synch_gen, 
			 rand_vertex, 
			 run_kla, 
			 run_dc, 
			 run_ds, 
			 run_cc, 
			 run_chaotic, 
			 run_ss_mis,
			 run_bucket_mis,
			 run_luby_mis,
			 k_level,
			 luby_algo,
			 cc_algo,
			 bucket_choice);
	}
	if(routing_choice == rt_hypercube) {
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_choice), reduction_choice, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_choice, amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  select_msg_gen(msg_gen, priority_msg_gen, trans, barrier_trans, g, distance, weight, thread_choice, flush_choice, eager_choice, delta_choice, dc_ds_choice, synch_gen, 
			 rand_vertex, 
			 run_kla, 
			 run_dc, 
			 run_ds, 
			 run_cc, 
			 run_chaotic, 
			 run_ss_mis,
			 run_bucket_mis,
			 run_luby_mis,
			 k_level,
			 luby_algo,
			 cc_algo,
			 bucket_choice);
	}
	if(routing_choice == rt_rook) {
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_choice), reduction_choice, amplusplus::rook_routing(trans.rank(), trans.size()));
	  MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), reduction_choice, amplusplus::rook_routing(trans.rank(), trans.size()));
	  select_msg_gen(msg_gen, priority_msg_gen, trans, barrier_trans, g, distance, weight, thread_choice, flush_choice, eager_choice, delta_choice, dc_ds_choice, synch_gen, 
			 rand_vertex, 
			 run_kla, 
			 run_dc, 
			 run_ds, 
			 run_cc, 
			 run_chaotic, 
			 run_ss_mis,
			 run_bucket_mis,
			 run_luby_mis,
			 k_level,
			 luby_algo,
			 cc_algo,
			 bucket_choice);
	}
      }
    }
 
    if(no_reductions) {
      if (trans.rank() == 0) std::cout << " Reduction: 0 " << "coalescing " << coalescing_choice << " pri coalescing " 
				       << priority_coalescing_size 
				       << " routing " << routing_choice
				       << " run_cc " << run_cc
				       << std::endl;
      if(routing_choice == rt_none) {
	typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(coalescing_choice)));
	MessageGenerator priority_msg_gen((amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1)));
	select_msg_gen(msg_gen, priority_msg_gen, trans, barrier_trans, g, distance, weight, thread_choice, flush_choice, eager_choice, delta_choice, dc_ds_choice, synch_gen, 
		       rand_vertex, 
		       run_kla, 
		       run_dc, 
		       run_ds, 
		       run_cc, 
		       run_chaotic, 
		       run_ss_mis,
		       run_bucket_mis,
		       run_luby_mis,
		       k_level,
		       luby_algo,
		       cc_algo,
		       bucket_choice);
      }
      if(routing_choice == rt_hypercube) {
	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_choice), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), amplusplus::hypercube_routing(trans.rank(), trans.size()));
	select_msg_gen(msg_gen, priority_msg_gen, trans, barrier_trans, g, distance, weight, thread_choice, flush_choice, eager_choice, delta_choice, dc_ds_choice, synch_gen, 
		       rand_vertex, 
		       run_kla, 
		       run_dc, 
		       run_ds, 
		       run_cc, 
		       run_chaotic, 
		       run_ss_mis,
		       run_bucket_mis,
		       run_luby_mis,
		       k_level,
		       luby_algo,
		       cc_algo,
		       bucket_choice);
      }
      if(routing_choice == rt_rook) {
	typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(coalescing_choice), amplusplus::rook_routing(trans.rank(), trans.size()));
	MessageGenerator priority_msg_gen(amplusplus::counter_coalesced_message_type_gen(priority_coalescing_size, 1), amplusplus::rook_routing(trans.rank(), trans.size()));
	select_msg_gen(msg_gen, priority_msg_gen, trans, barrier_trans, g, distance, weight, thread_choice, flush_choice, eager_choice, delta_choice, dc_ds_choice, synch_gen, 
		       rand_vertex, 
		       run_kla, 
		       run_dc, 
		       run_ds, 
		       run_cc, 
		       run_chaotic, 
		       run_ss_mis,
		       run_bucket_mis,
		       run_luby_mis,
		       k_level,
		       luby_algo,
		       cc_algo,
		       bucket_choice);
      }
    }
  }
  
  template<typename MessageGenerator, typename Digraph, typename DistanceMap, typename WeightMap, typename RNG, typename RandVertex>
  void select_msg_gen(MessageGenerator &msg_gen, 
		      MessageGenerator &priority_msg_gen, 
		      amplusplus::transport &trans, 
		      amplusplus::transport &barrier_trans, 
		      Digraph &g, 
		      DistanceMap &distance, 
		      WeightMap &weight, 
		      unsigned int num_threads, 
		      unsigned int flushFreq, 
		      unsigned int eager_limit, 
		      unsigned int delta, 
		      std::string& dc_ds, 
		      RNG synch_gen, 
		      RandVertex rand_vertex, 
		      bool run_kla, 
		      bool run_dc, 
		      bool run_ds, 
		      bool run_cc, 
		      bool run_chaotic, 
		      bool run_ss_mis,
		      bool run_bucket_mis,
		      bool run_luby_mis,
		      size_t k_level,
		      std::string luby_algo,
		      std::string cc_algo,
		      size_t bucket_choice) {
    //run_ds = !run_dc;

    typedef typename graph_traits<Digraph>::vertex_descriptor Vertex;

    time_type total_time = 0;
    std::vector<time_type> all_times(num_sources);

    std::vector<state_t> misvec(num_vertices(g), MIS_UNFIX);
    typedef typename property_map<Digraph, vertex_index_t>::type VertexIndexMap;
    typedef iterator_property_map<typename std::vector<state_t>::iterator, VertexIndexMap>  MISMap;
    MISMap mis(misvec.begin(), get(vertex_index, g));


#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    clear_cumulative_buffer_stats();
#endif
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    cumulative_hits = 0; cumulative_tests = 0;
#endif
#ifdef PBGL2_PRINT_WORK_STATS
    clear_cumulative_work_stats();
#endif


    // For GTEPS
    double total_elapsed_time = 0.0;
    size_t *edge_traversed =(size_t *) calloc(num_sources, sizeof(size_t));
    double *elapsed_time = (double *) calloc(num_sources, sizeof(double));

    for (unsigned int source_i = 0; source_i < num_sources; ++source_i) {

#ifdef PRINT_ET
      if (trans.rank() == 0)
        std::cout << print_time(get_time() - job_start) << "s elapsed since job start\n";
#endif
      Vertex current_source = vertex(rand_vertex(synch_gen), g);
      time_type time;
      
      if(run_dc) {
	if(dc_ds == "vector") {
	  // Note : We need to modify vector data srtucture to accomadate new interface
	  std::cout << " This (vector data structure) option is not supported." 
		    << " We need to change vector implementation to satisfy new interface"
		    << std::endl;
	  exit(1);
	  /*
	  time = run_distributed_control<boost::graph::distributed::vector_of_vector_gen>
	    (trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq, eager_limit);
	  */
	} else if (dc_ds == "nodeq") { // priority queue per node
	  time = run_distributed_control_node<boost::graph::distributed::node_priority_queue_gen>
	    (trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq, eager_limit);
	} else if (dc_ds == "ms_nodeq"){ // priority queue per node
	  time = run_distributed_control_node<boost::graph::distributed::ms_node_priority_queue_gen>
	    (trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq, eager_limit);
	} else if (dc_ds == "numaq") { // priority queue per numa node
	  time = run_distributed_control_node<boost::graph::distributed::numa_priority_queue_gen>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit,
	     true); // numa true
	} else if (dc_ds == "ms_numaq") { // priority queue per numa node
	  time = run_distributed_control_node<boost::graph::distributed::ms_numa_priority_queue_gen>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit,
	     true); // numa true
	} else if (dc_ds == "pheet64q"){ 
	  time = run_distributed_control_pheet<boost::graph::distributed::pheet_priority_queue_gen_64>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit); 
	}  else if (dc_ds == "pheet128q"){ 
	  time = run_distributed_control_pheet<boost::graph::distributed::pheet_priority_queue_gen_128>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit); 
	}  else if (dc_ds == "pheet256q"){ 
	  time = run_distributed_control_pheet<boost::graph::distributed::pheet_priority_queue_gen_256>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit); 
	}  else if (dc_ds == "pheet512q"){ 
	  time = run_distributed_control_pheet<boost::graph::distributed::pheet_priority_queue_gen_512>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit); 
	}  else if (dc_ds == "pheet1024q"){ 
	  time = run_distributed_control_pheet<boost::graph::distributed::pheet_priority_queue_gen_1024>
	    (trans, barrier_trans, g, weight, distance, 
	     current_source, num_threads, n, 
	     verify, msg_gen, priority_msg_gen, 
	     flushFreq, eager_limit); 
	} else if (dc_ds == "threadq"){
	  time = run_distributed_control(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, priority_msg_gen, flushFreq, eager_limit);
	} else if (dc_ds == "buffer") {
	  std::cout << "Ignoring buffer data structure for Distirbuted Cotrol ..." << std::endl;
	} else {
	  std::cout << "Data structure not specified ! " << std::endl;
	  assert(false);
	}
      } else if (run_cc) { 

	if (cc_algo == "run_dd_cc") {
	  time = run_data_driven_cc(trans, 
				    barrier_trans, 
				    g, 
				    num_threads, n, 
				    verify, 
				    print_cc_stats,
				    allstats,
				    msg_gen, 
				    priority_msg_gen, flushFreq, eager_limit);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_sv_cc") {
	  time = run_sv_c_c(trans, 
			    g, 
			    num_threads, 
			    verify, 
			    false, // level_sync 
			    msg_gen,
			    vertices_per_lock);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_sv_cc_level_sync") {
	  time = run_sv_c_c(trans, 
			    g, 
			    num_threads, 
			    verify, 
			    true, // level_sync 
			    msg_gen,
			    vertices_per_lock);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_sv_cc_opt") {
	  time = run_sv_cc_optimized(trans, 
			    g, 
			    num_threads, 
			    verify, 
			    false, // level_sync 
			    msg_gen,
			    vertices_per_lock);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_sv_cc_opt_level_sync") {
	  time = run_sv_cc_optimized(trans, 
			    g, 
			    num_threads, 
			    verify, 
			    true, // level_sync 
			    msg_gen,
			    vertices_per_lock);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}


	if (cc_algo == "run_sv_ps_cc") {
	  time = run_ps_sv_cc(trans, 
			      g, 
			      num_threads, 
			      verify, 
			      false, // level_sync 
			      msg_gen,
			      vertices_per_lock);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_sv_ps_cc_level_sync") {
	  time = run_ps_sv_cc(trans, 
			      g, 
			      num_threads, 
			      verify, 
			      true, // level_sync 
			      msg_gen,
			      vertices_per_lock);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_cc_chaotic") {
	  time = run_chaotic_cc(trans,
				barrier_trans,
				g, 
				num_threads, 
				n, 
				verify, 
				print_cc_stats,
				allstats,
				msg_gen);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_cc_ds") {
	  if (id_distribution == vertical)
	    time = run_delta_cc(trans,
				barrier_trans,
				g, 
				num_threads, 
				n, 
				bucket_choice,
				block_id_distribution<Digraph>(g, n),
				verify, 
				print_cc_stats,
				allstats,
				msg_gen);
	  else
	    time = run_delta_cc(trans,
				barrier_trans,
				g, 
				num_threads, 
				n, 
				bucket_choice,
				row_id_distribution<Digraph>(g, trans.size()),
				verify, 
				print_cc_stats,
				allstats,
				msg_gen);


	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}

	if (cc_algo == "run_level_sync_cc") {
	  time = run_level_sync_cc(trans,
				   barrier_trans,
				   g, 
				   num_threads, 
				   n, 
				   verify, 
				   print_cc_stats,
				   allstats,
				   msg_gen);

	  { amplusplus::scoped_epoch epoch(barrier_trans); }
	}



      } else if(run_ds) { // run_ds
	// For now we don't have the parameter for level_sync.  It can be added later on.
	if (dc_ds == "buffer") {
	  time = run_delta_stepping(trans, 
				    barrier_trans, 
				    g, 
				    weight, 
				    distance, 
				    current_source, 
				    delta, 
				    num_threads, n, verify, false, msg_gen);
	} else if (dc_ds == "nodeq") {
	  time = run_delta_stepping_node<boost::graph::distributed::node_priority_queue_gen>(trans, 
											     barrier_trans, g, 
											     weight, distance, 
											     current_source, delta, 
											     num_threads, n, verify, false, msg_gen,
											     flushFreq);
	} else if (dc_ds == "ms_nodeq") {
	  time = run_delta_stepping_node<boost::graph::distributed::ms_node_priority_queue_gen>(trans, 
												barrier_trans, g, 
												weight, distance, 
												current_source, delta, 
												num_threads, n, verify, false, msg_gen,
												flushFreq);
	} else if (dc_ds == "threadq") {
	  time = run_delta_stepping_thread<boost::graph::distributed::thread_priority_queue_gen>(trans, 
												 barrier_trans, g, 
												 weight, distance, 
												 current_source, delta, 
												 num_threads, n, verify, false, msg_gen, 
												 flushFreq);
	} else if (dc_ds == "mtthreadq") {
	  time = run_delta_stepping_thread<boost::graph::distributed::multi_priority_queue_gen>(trans, 
												 barrier_trans, g, 
												 weight, distance, 
												 current_source, delta, 
												 num_threads, n, verify, false, msg_gen, 
												 flushFreq);
	} else if (dc_ds == "numaq") {
	  time = run_delta_stepping_numa<boost::graph::distributed::numa_priority_queue_gen>(trans, 
											     barrier_trans, g, 
											     weight, distance, 
											     current_source, delta, 
											     num_threads, n, verify, false, msg_gen,
											     flushFreq);
	} else if (dc_ds == "ms_numaq") {
	  time = run_delta_stepping_numa<boost::graph::distributed::ms_numa_priority_queue_gen>(trans, 
												barrier_trans, g, 
												weight, distance, 
												current_source, delta, 
												num_threads, n, verify, false, msg_gen,
												flushFreq);
	} else {
	  std::cout << "Data structure not specified ! " << std::endl;
	  assert(false);
	}
      } else if(run_chaotic) {
	// chaotic SSSP
	time = run_distributed_control_chaotic(trans, barrier_trans, g, weight, distance, current_source, num_threads, n, verify, msg_gen, flushFreq, eager_limit);
      } else if(run_kla) {
	if (dc_ds == "buffer") {
	  time = run_kla_sssp(trans, barrier_trans, 
			      g, weight, distance, 
			      current_source, 
			      num_threads, n, verify, 
			      false, msg_gen, k_level);
	} else if (dc_ds == "numaq") {
	  time = run_kla_sssp_numa<boost::graph::distributed::numa_priority_queue_gen>(trans, 
										       barrier_trans, 
										       g, weight, distance, 
										       current_source, 
										       num_threads, n, verify, 
										       false, msg_gen, k_level, flushFreq);
	} else if (dc_ds == "ms_numaq") {
	  time = run_kla_sssp_numa<boost::graph::distributed::ms_numa_priority_queue_gen>(trans, 
											  barrier_trans, 
											  g, weight, distance, 
											  current_source, 
											  num_threads, n, verify, 
											  false, msg_gen, k_level,
											  flushFreq);
	} else if (dc_ds == "nodeq") {
	  time = run_kla_sssp_node<boost::graph::distributed::node_priority_queue_gen>(trans, 
										       barrier_trans, 
										       g, weight, distance, 
										       current_source, 
										       num_threads, n, verify, 
										       false, msg_gen, k_level, flushFreq);
	} else if (dc_ds == "ms_nodeq") {
	  time = run_kla_sssp_node<boost::graph::distributed::ms_node_priority_queue_gen>(trans, 
											  barrier_trans, 
											  g, weight, distance, 
											  current_source, 
											  num_threads, n, verify, 
											  false, msg_gen, k_level, flushFreq);
	} else if (dc_ds == "threadq"){
	  time = run_kla_sssp_thread<boost::graph::distributed::thread_priority_queue_gen>(trans, 
											   barrier_trans, 
											   g, weight, distance, 
											   current_source, 
											   num_threads, n, verify, 
											   false, msg_gen, k_level,
											   flushFreq);
	} else {
	  assert(false);
	}
      } else if(run_ss_mis) {
	time = run_fix_mis<boost::graph::distributed::thread_priority_queue_gen>(trans,
										 barrier_trans,
										 g,
										 mis,
										 num_threads,
										 n,
										 verify,
										 msg_gen,
										 flushFreq);
      } else if(run_bucket_mis) {
	time = run_fix_mis_bucket(trans,
				 barrier_trans,
				 g,
				 mis,
				 num_threads,
				 n,
				 verify,
				 msg_gen,
				 flushFreq);
      } else if(run_luby_mis) {
	if (luby_algo == "A") {
	  time = run_luby_maximal_is<boost::graph::distributed::select_a_functor_gen>(trans,
										      barrier_trans,
										      g,
										      mis,
										      num_threads,
										      n,
										      verify,
										      msg_gen,
										      flushFreq);

	} else if (luby_algo == "AV1") {
	  time = run_luby_maximal_is<boost::graph::distributed::select_a_vertex_functor_gen>(trans,
											     barrier_trans,
											     g,
											     mis,
											     num_threads,
											     n,
											     verify,
											     msg_gen,
											     flushFreq);

	} else if (luby_algo == "AV2"){
	  time = run_luby_maximal_is<boost::graph::distributed::select_a_v2_functor_gen>(trans,
											     barrier_trans,
											     g,
											     mis,
											     num_threads,
											     n,
											     verify,
											     msg_gen,
											     flushFreq);

	} else if (luby_algo == "B") {
	  time = run_luby_maximal_is<boost::graph::distributed::select_b_functor_gen>(trans,
										      barrier_trans,
										      g,
										      mis,
										      num_threads,
										      n,
										      verify,
										      msg_gen,
										      flushFreq);

	} else {
	  std::cerr << "Invalid algorithm type for Luby MIS."
		    << " Available algorithms are A, AV1, AV2 and B"
		    << std::endl;
	  assert(false);
	}


      } else {
	assert(false);
      }

#ifdef CALCULATE_GTEPS
      if (!run_cc) {
	elapsed_time[source_i] = time;
	edge_traversed[source_i] = get_gteps(trans, g, distance, weight);
      }
#endif
      if (time == -1.) { // Not enough vertices visited
	--source_i; continue;
      }

      total_time += time;
      all_times[source_i] = time;

    }

    if (trans.rank() == 0) {
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
      unsigned long long buffers{}, messages{};
      for(unsigned int i = 0; i < cumulative_flushes.size(); ++i) {
	buffers += cumulative_flushes[i];
	messages += cumulative_flushes[i] * i;
      }

      float msg_per_buf = -1;
      if (buffers != 0) {
	msg_per_buf = messages/buffers;
      }
#endif

      if(run_dc) {
	std::cout << "Total Distributed Control (DC)-(" << dc_ds 
		  << ") time for " << num_sources << " sources = "
		  << print_time(total_time) << " ("
		  << print_time(total_time / num_sources) << " per source), "
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
		  << "messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), "
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
      } else if(run_chaotic) {
	std::cout << "Total Chaotic SSSP " 
		  << "time for " << num_sources << " sources = "
		  << print_time(total_time) << " ("
		  << print_time(total_time / num_sources) << " per source), "
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
		  << "messages = " << cumulative_messages << " (" << (cumulative_messages / num_sources) << " per source), "
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
      } else if(run_ds) { // run_ds
	std::cout << "Total Delta-Stepping (DS-" << dc_ds
		  << ") time for " << num_sources << " sources = " 
		  << print_time(total_time) << " (" 
		  << print_time(total_time / num_sources) << " per source), delta = " << delta 
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
	  PBGL2_DS_PRINT // The print macro for work stats
#endif
		  << std::endl;
      } /*else if (run_cc) { // run connected components
#ifdef CC_ALGORITHMS
	std::cout << "Total data driven cc time for " << num_sources << " sources = " 
		  << print_time(total_cc_dd_time) << " (" 
		  << print_time(total_cc_dd_time / num_sources) << " per source)"
		  << " Total shiloach-vishkin time " << print_time(total_cc_sv_time) << " ("
		  << print_time(total_cc_sv_time / num_sources) << " per source)"
		  << " Total parallel shiloach-vishkin time " << print_time(total_cc_sv_ps_time) << " ("
		  << print_time(total_cc_sv_ps_time / num_sources) << " per source)"
		  << " Total local difference time " << print_time(total_cc_ld_time) << " ("
		  << print_time(total_cc_ld_time / num_sources) << " per source)."
		  << " Running Scale : " << scale
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
	  PBGL2_DS_PRINT // The print macro for work stats
#endif
		  << std::endl;
#endif // end of CC_ALGORITHMS
}*/ else if(run_kla) { // run_kla
	std::cout << "Total kla-sssp time for " << num_sources << " sources = " 
		  << print_time(total_time) << " (" 
		  << print_time(total_time / num_sources) << " per source), klevel = " 
		  << k_level << "%"
#ifdef PBGL2_PRINT_WORK_STATS
	  PBGL2_DS_PRINT // The print macro for work stats
#endif
		  << std::endl; 
      }

      time_type mean = total_time / num_sources;
      time_type min = 0;
      time_type q1 = 0;
      time_type median = 0;
      time_type q3 = 0;
      time_type max = 0;
      time_type stddev = 0;

      time_statistics(all_times, 
		      mean, 
		      min,
		      q1,
		      median,
		      q3,
		      max,
		      stddev);

      std::cout << "MEAN : " << print_time(mean)
		<< " STDDEV : " << print_time(stddev)
		<< " MIN : " << print_time(min)
		<< " Q1 : " << print_time(q1)
		<< " MEDIAN : " << print_time(median)
		<< " Q3 : " << print_time(q3)
		<< " MAX : " << print_time(max)
		<< std::endl;

#ifdef CALCULATE_GTEPS
      std::cout << "==============TEPS Statistics=================" << std::endl; 
      double *tm = (double*)malloc(sizeof(double)*num_sources);
      double *stats = (double*)malloc(sizeof(double)*9);

      // TODO We need to fix stat when num_sources < 4
      if (!run_cc && num_sources > 4) {
	for(int i = 0; i < num_sources; i++)
	  tm[i] = edge_traversed[i]/elapsed_time[i];

	statistics (stats, tm, num_sources);
	PRINT_GRAPH500_STATS("TEPS", 1);
      }

      free(tm);
      free(stats);
#endif
      free(edge_traversed);
      free(elapsed_time);
    }
  }

  template<typename Graph500Iter, typename Routing, typename Distribution, typename Trans, typename EdgeSize, typename UInt, typename Edges, typename WeightIterator>
  void do_distribute(Distribution &distrib, Trans &trans, EdgeSize e_start, EdgeSize e_count, UInt a, UInt b, Edges &edges, WeightIterator weight_iter,
		     int num_threads) {
    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, Routing> MessageGenerator;
    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(distribution_coalescing_size), Routing(trans.rank(), trans.size()));

    time_type start = get_time();

    // Multithreaded distribution
    trans.set_nthreads(num_threads);

#ifdef PRINT_DEBUG
    std::cout << "Thread start : " << e_start << "Total edge count : " << e_count
	      << " for rank :" << trans.rank()
	      << std::endl;
#endif

    Edges* thread_edges = new Edges[num_threads];
    EdgeSize std_count_per_thread = (e_count + num_threads - 1) / num_threads;
    //   EdgeSize remainder = e_count % num_threads;

#ifdef PRINT_DEBUG
    std::cout << "std_count_per_thread : " << std_count_per_thread << std::endl;
#endif

    EdgeSize thread_start = e_start;
    EdgeSize count_per_thread = std_count_per_thread;

    // create message type
    typedef typename std::iterator_traits<Graph500Iter>::value_type value_type;
    typedef detail::vector_write_handler<Edges> vec_write_handler;

    // not sure we need this here
    amplusplus::register_mpi_datatype<value_type>();
    
    typedef typename MessageGenerator::template call_result<value_type, vec_write_handler,
      detail::distrib_to_pair_owner_map_t<Distribution>, amplusplus::no_reduction_t>::type 
      write_msg_type;

    write_msg_type write_msg(msg_gen, trans, detail::distrib_to_pair_owner_map<Distribution>(distrib), 
			     amplusplus::no_reduction);

    // set the handler (must be in main thread)
    write_msg.set_handler(vec_write_handler(thread_edges));
    threaded_distributor<write_msg_type> td(write_msg);

    boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
    for (int i = 0; i < num_threads - 1; ++i) {
#ifdef PRINT_DEBUG
      if (_RANK == 0) {
	std::cout << "Thread : " << i+1 << " start : " << thread_start
		  << " thread end : " << (thread_start + count_per_thread)
		  << std::endl;
      }
#endif
      boost::thread thr(boost::ref(td), i+1, Graph500Iter(scale, thread_start, a, b, weight_iter),
			Graph500Iter(scale, thread_start + count_per_thread, a, b, weight_iter),
			flip_pair(), trans);
      thread_start = thread_start + count_per_thread;
      threads[i].swap(thr);
    }

    // Allocate remainder to main thread
    count_per_thread = std::min(std_count_per_thread, (e_count-thread_start));

#ifdef PRINT_DEBUG
    std::cout << "Thread : " << 0 << " start : " << thread_start
		<< " thread end : " << (thread_start + count_per_thread)
		<< std::endl;
#endif

    td(0, Graph500Iter(scale, thread_start, a, b, weight_iter),
       Graph500Iter(scale, thread_start + count_per_thread, a, b, weight_iter),
       flip_pair(), trans);

    // Wait till all threads finish
    for (int i = 0; i < num_threads - 1; ++i)
      threads[i].join();

    // Merge all edges - in threads
    //    time_type m_start = get_time();

    typename Edges::size_type total_sz = 0;
    std::vector<typename Edges::size_type> positions(num_threads);
    for (int i = 0; i < num_threads; ++i) {
      positions[i] = total_sz;
      total_sz += thread_edges[i].size();
    }

    edges.resize(total_sz);
    threaded_merger<Edges> tm(edges);

    boost::scoped_array<boost::thread> m_threads(new boost::thread[num_threads - 1]);
    for (int i = 0; i < num_threads - 1; ++i) {
      typename Edges::iterator pos_ite = edges.begin();
      std::advance(pos_ite, positions[i+1]);
      boost::thread thr(boost::ref(tm), pos_ite, 
			thread_edges[i+1].begin(), thread_edges[i+1].end());
      m_threads[i].swap(thr);
    }

    tm(edges.begin(), thread_edges[0].begin(), thread_edges[0].end());

    // Wait till all threads finish
    for (int i = 0; i < num_threads - 1; ++i)
      m_threads[i].join();

    //time_type m_end = get_time();

    trans.set_nthreads(1);
    delete[] thread_edges;
    { amplusplus::scoped_epoch epoch(trans); }
    time_type end = get_time();
    std::cout << "Graph distribution permutation took " << (end-start)
	      << " time" << std::endl;

    // sort edges in parallel
    //__gnu_parallel::sort(edges.begin(), edges.end());
  }
};

void q_test() {
  typedef std::pair<int, int> vertex_distance_data;
  struct default_comparer {
    bool operator()(const vertex_distance_data& vd1, const vertex_distance_data& vd2) {
      return vd1.second > vd2.second;
    }
  };

  cds::Initialize();
  {
    // Initialize Hazard Pointer singleton
    cds::gc::HP hpGC(72);
    // Attach for the main thread
    cds::threading::Manager::attachThread();

    //  template<typename vertex_distance, typename Compare>
    /*  struct reverse_compare {
	private:
	default_comparer c;
	public:
	inline bool operator()(const vertex_distance_data& vd1, const vertex_distance_data& vd2) {
	return c(vd2, vd1);
	}
	};*/

    graph::distributed::ms_node_priority_queue<vertex_distance_data, default_comparer> npq(10);
    graph::distributed::node_priority_queue<vertex_distance_data, default_comparer> pq(10);
    //    graph::distributed::ellen_bin_priority_queue<vertex_distance_data, default_comparer> epq(10);
    std::priority_queue<vertex_distance_data, std::vector<vertex_distance_data>, default_comparer> dpq;
    vertex_distance_data vd1(10, 22);
    vertex_distance_data vd2(1, 23);
    vertex_distance_data vd3(5, 34);
    vertex_distance_data vd4(8, 12);
    vertex_distance_data vd5(11, 210);
    vertex_distance_data vd6(7, 2);
    vertex_distance_data vd7(8, 650);
    vertex_distance_data vd8(5, 30);
    vertex_distance_data vd9(2, 22);
    vertex_distance_data vd10(9, 21);

    npq.put(vd1,0);
    npq.put(vd2,1);
    npq.put(vd3,0);
    npq.put(vd4,2);
    npq.put(vd5,3);
    npq.put(vd6,4);
    npq.put(vd7,5);
    npq.put(vd8,6);
    npq.put(vd9,2);
    npq.put(vd10,4);


    pq.put(vd1,0);
    pq.put(vd2,1);
    pq.put(vd3,0);
    pq.put(vd4,2);
    pq.put(vd5,3);
    pq.put(vd6,4);
    pq.put(vd7,5);
    pq.put(vd8,6);
    pq.put(vd9,2);
    pq.put(vd10,4);


    /*    epq.put(vd1,0);
    epq.put(vd2,1);
    epq.put(vd3,0);
    epq.put(vd4,2);
    epq.put(vd5,3);
    epq.put(vd6,4);
    epq.put(vd7,5);
    epq.put(vd8,6);
    epq.put(vd9,2);
    epq.put(vd10,4);*/


    dpq.push(vd1);
    dpq.push(vd2);
    dpq.push(vd3);
    dpq.push(vd4);
    dpq.push(vd5);
    dpq.push(vd6);
    dpq.push(vd7);
    dpq.push(vd8);
    dpq.push(vd9);
    dpq.push(vd10);

    std::cout << "========== cds :: mspq ===========" << std::endl;
    // should give the shortest distance first
    vertex_distance_data vdx1;
    while(npq.pop(vdx1, 0)) {
      std::cout << vdx1.second << std::endl;
    }
    assert(npq.empty());
  
    std::cout << "========== cds :: fcpq ===========" << std::endl;
    // should give the shortest distance first
    vertex_distance_data vdx2;
    while(pq.pop(vdx2, 0)) {
      std::cout << vdx2.second << std::endl;
    }
    assert(pq.empty());

    /*    std::cout << "========== cds :: ellenbin ===========" << std::endl;
    // should give the shortest distance first
    vertex_distance_data vdx3;
    while(epq.pop(vdx3, 0)) {
      std::cout << vdx3.second << std::endl;
    }
    assert(epq.empty());*/

    std::cout << "========== std::queue===========" << std::endl;
    while(!dpq.empty()) {
      std::cout << dpq.top().second << std::endl;
      dpq.pop();
    }

    cds::threading::Manager::detachThread();
  }
  cds::Terminate();
  exit(0);

}

int main(int argc, char* argv[]) {
  
  //  q_test();
  
  dc_test t(argc, argv);
  t.test_main(argc, argv);
}


