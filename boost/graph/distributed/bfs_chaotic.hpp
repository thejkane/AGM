// Copyright 2017 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_GRAPH_BFS_CHAOTIC
#define BOOST_GRAPH_BFS_CHAOTIC

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
#include <boost/graph/distributed/owner_defs.hpp>
#include <queue>
#include <algorithm> // for std::min, std::max
#include <iostream>
#include <atomic>
#include <tuple>
#include <climits>
#include <time.h>

typedef uint32_t Distance;
// work stats - <useful, invalidated, rejected>
typedef std::tuple<uint64_t, uint64_t, uint64_t> work_stats_t;

namespace boost {
namespace graph {
namespace distributed {

template<typename Graph,
	 typename DistanceMap,
	 typename MessageGenerator = amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class bfs_chaotic {

  // Type definitions for Vertex and Distance
  // Vertex information is in graph_traits::vertex_descriptor
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  // This is going to be our message
  typedef std::pair<Vertex, std::pair<std::pair<Vertex, Distance >, int64_t> > work_item_t;

  // Owner map type - To calculate which vertex belong to which node
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  // The handler, forward decleration
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

  struct empty_deleter {
  private:
    work_item_t *pwi;
  public:
    empty_deleter(work_item_t* p) : pwi(p){}
    typedef void result_type;
    void operator()() const {
      delete pwi;
    }
    //template <typename T> void operator()(const T&) const {}
  };



#ifndef NO_COALESCING
  typedef typename MessageGenerator::template call_result<
    work_item_t,
    processing_function,
    owner_from_pair<OwnerMap, work_item_t>,
    amplusplus::idempotent_combination_t<minimum_pair_first> >::type
  RelaxMessage;
#endif

public:

  bfs_chaotic(Graph& g,
	      DistanceMap distance_map,
	      amplusplus::transport &t,
	      int offset,
	      MessageGenerator message_gen =
	      MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
    graph(g),
    transport(t),
    core_offset(offset),
    distance_map(distance_map),
    owner(get(vertex_owner, g)),
#ifdef NO_COALESCING
    mt(transport.create_message_type<work_item_t>(0))
#else
    relax_msg(message_gen,
	      transport,
	      owner_from_pair<OwnerMap, work_item_t>(owner),
	      amplusplus::idempotent_combination(minimum_pair_first()))
#endif
  {
    initialize();
  }

#ifdef COLLECT_STATS
  void print_stats() {
    uint64_t totuseful = 0;
    uint64_t totinvalidated = 0;
    uint64_t totrejected = 0;

    for(unsigned int i = 0; i < useful.size(); ++i) {
      totuseful += useful[i];
      totinvalidated += invalidated[i];
      totrejected += rejected[i];
    }

    uint64_t g_totuseful = 0;
    uint64_t g_totinvalidated = 0;
    uint64_t g_totrejected = 0;

    MPI_Reduce(&totuseful, &g_totuseful, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totinvalidated, &g_totinvalidated, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totrejected, &g_totrejected, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);


    if (_RANK == 0) {
      std::cout << "Useful : " << g_totuseful << std::endl;
      std::cout << "Invalidated : " << g_totinvalidated << std::endl;
      std::cout << "Rejected : " << g_totrejected << std::endl;
    }
  }

  work_stats_t get_local_work_stats() {
    work_stats_t stats{ 0ul, 0ul, 0ul };
    uint64_t totuseful = 0;
    for(unsigned int i = 0; i < useful.size(); ++i) {
      totuseful += useful[i];
      std::get<1>(stats) += invalidated[i];
      std::get<2>(stats) += rejected[i];
    }

    std::get<0>(stats) = (totuseful - std::get<1>(stats));
    return stats;
  }
#endif

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  typedef std::pair<unsigned long long, unsigned long long> cache_stats;
  cache_stats
  get_cache_stats() {
    cache_stats stats(0, 0);
#ifndef NO_COALESCING
    for(size_t i = 0; i < relax_msg.counters.hits.size(); ++i) {
      stats.first += relax_msg.counters.hits[i];
      stats.second += relax_msg.counters.tests[i];
    }
#endif
    return stats;
  }
#endif // AMPLUSPLUS_PRINT_HIT_RATES

#ifdef PBGL2_PRINT_WORK_STATS
  // useful, invalidated, useless, rejected
  typedef std::tuple<unsigned long long, unsigned long long, unsigned long long, unsigned long long> work_stats_t;
  work_stats_t
  get_work_stats() {
    work_stats_t stats{ 0ul, 0ul, 0ul, 0ul };
    for(unsigned int i = 0; i < useful.size(); ++i) {
      std::get<0>(stats) += useful[i];
      std::get<1>(stats) += invalidated[i];
      std::get<2>(stats) += useless[i];
      std::get<3>(stats) += rejected[i];
    }
    return stats;
  }
#endif // PBGL2_PRINT_WORK_STATS

  void set_source(Vertex s) { source = s; }
  void set_abort_level(int level) { abort_level = level; }
  // What is actually getting called is this operator
  void operator() (int tid) { run(source, tid); }
  void run(Vertex s, int tid = 0);
  time_type get_start_time() { return start_time; }

  uint64_t get_traversed_edges() {
    uint64_t edges = 0;
    for (int i=0; i < edges_traversed.size(); ++i) {
      edges += edges_traversed[i];
    }

    return edges;
  }

private:

  int64_t get_current_time() {
    auto now = std::chrono::system_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::nanoseconds>(now);
    auto value = now_ms.time_since_epoch();
    return value.count();
  }

  int64_t get_duration_ms(int64_t t) {
    auto current = get_current_time();
    return (current - t);
  }

  void initialize();
  void relax(work_item_t, int tid);
  void handle_queue(const int, amplusplus::transport::end_epoch_request&);

  const int dummy_first_member_for_init_order; // Unused
  const Graph& graph;
  amplusplus::transport& transport;
  int core_offset;
  DistanceMap distance_map;
  const OwnerMap owner;
  Vertex source;
#ifdef NO_COALESCING
  typedef amplusplus::message_type<work_item_t> mt_t;
  mt_t mt;
  //std::vector<mt_t*> msg_types;
#else
  RelaxMessage relax_msg;
#endif
  shared_ptr<amplusplus::detail::barrier> t_bar;
  unsigned int in_edge_factor;
  time_type start_time;
  std::vector<uint64_t> edges_traversed;
#ifdef COLLECT_STATS
  std::vector<unsigned long long> useful, invalidated, rejected;
  bool first_relax;
#endif // PBGL2_PRINT_WORK_STATS
  int abort_level;
};

#define BFS_CHAOTIC_PARMS					\
  typename Graph, typename DistanceMap, typename MessageGenerator

#define BFS_CHAOTIC_TYPE						\
  bfs_chaotic<Graph, DistanceMap, MessageGenerator>

template<BFS_CHAOTIC_PARMS> void
BFS_CHAOTIC_TYPE::initialize() {
  transport.template downcast_to_impl<amplusplus::mpi_transport_event_driven>()->set_use_ssend(true);

#ifdef COLLECT_STATS
  useful.resize(transport.get_nthreads(), 0);
  invalidated.resize(transport.get_nthreads(), 0);
  rejected.resize(transport.get_nthreads(), 0);
  first_relax = true;
#endif // COLLECT_STATS

  edges_traversed.resize(transport.get_nthreads(), 0);

#ifdef NO_COALESCING
  /*  msg_types.resize(transport.get_nthreads());
  for (int i=0; i < transport.get_nthreads(); ++i) {
    msg_types[i] = new mt_t(transport.create_message_type<work_item_t>(0));
    msg_types[i]->set_max_count(1);
    msg_types[i]->set_handler(processing_function(*this));
    }*/


  //amplusplus::valid_rank_set possible_dests_ = amplusplus::valid_rank_set();
  //amplusplus::valid_rank_set possible_sources_ = amplusplus::valid_rank_set();
  mt.set_max_count(1);
  mt.set_handler(processing_function(*this));
  //mt.set_possible_dests(possible_dests_);
  //mt.set_possible_sources(possible_sources_);
#else
  relax_msg.set_handler(processing_function(*this));
#endif

}

template<BFS_CHAOTIC_PARMS> void
BFS_CHAOTIC_TYPE::run(Vertex s, int tid) {

  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    if (0 == tid) {
      int nthreads = transport.get_nthreads();

      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    }

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

    t_bar->wait();

#ifdef CRAYPAT
    if (PAT_region_begin ( 1, "tcrun" ) == PAT_API_FAIL) {
      std::cout << "PAT begin failed ! " << std::endl;
      assert(false);
    }
#endif

    start_time = get_time();
    // Start the algorithm
    {
      amplusplus::scoped_epoch epoch(transport); 
      // If this is main thread and if node is the owner of the vertex do relax
      // for source with distance 0
      s = 2;
      if (get(owner, s) == transport.rank() && tid == 0) {
	std::cout << "source : " << s << std::endl;
	relax(work_item_t(s, std::make_pair(std::make_pair(s, 0), get_current_time())), tid);
      }
    }

  }
}

template<BFS_CHAOTIC_PARMS> void
BFS_CHAOTIC_TYPE::relax(work_item_t vd, int tid) {
  using boost::parallel::val_compare_and_swap;

  // handler is called per each edge
  (edges_traversed[tid])++;

#ifdef CYCLE_TIME
  auto valstart =  clock();
  auto valend =  clock(); 
  auto relstart =  clock();
  auto relend =  clock();
  auto ownstart =  clock();
  auto ownend =  clock();
  auto sendstart =  clock();
  auto sendend =  clock();
  auto start =  clock();
#else
  uint64_t totalsend = 0;
  unsigned int totalsendcnt = 0;
  auto valstart = std::chrono::high_resolution_clock::now();
  auto valend = std::chrono::high_resolution_clock::now();
  auto relstart = std::chrono::high_resolution_clock::now();
  auto relend = std::chrono::high_resolution_clock::now();
  auto ownstart = std::chrono::high_resolution_clock::now();
  auto ownend = std::chrono::high_resolution_clock::now();
  auto sendstart = std::chrono::high_resolution_clock::now();
  auto sendend = std::chrono::high_resolution_clock::now();
  auto start = std::chrono::high_resolution_clock::now();
#endif

  amplusplus::transport::rank_type last_owner;

  Vertex v = vd.first;
  Distance new_distance = vd.second.first.second;
  Vertex s = vd.second.first.first;


  std::string color = "\033["+std::to_string(0)+"m";
  if (new_distance > 0) {
    color = "\033[0;"+std::to_string(new_distance+29)+"m";
  }

  int64_t senttime = vd.second.second;
  std::clog << color << "[[" << tid << "](" << v << ", " << s << ")" << new_distance << "] " << std::flush; 
  //std::cout << "duration : " << get_duration_ms(senttime) << " ns" << std::endl;

#ifdef COLLECT_STATS
  if (!first_relax && (v == s)) {
    std::cout << "Self edges present !" << std::endl;
    assert(false);
  }

  if (first_relax)
    first_relax = false;
#endif

  //assert(v != s);

  // cut it out if distance is 2
  if (new_distance == abort_level)
    return;

  // we are processing an incoming edge for vertex v
  // increase incoming edge count
  //++processed_in_edge_count_map[v];

  // get existing distance
  // need last_old_dist to atomically update distance map
  auto readmapts = get_current_time();
  Distance old_dist = get(distance_map, v), last_old_dist;
  auto readmapte = get_current_time();

  //  std::cout << "read time : " << (readmapte - readmapts) << std::endl;


  // check whether we got a better distance if we got better
  // distance update the distance map
  // credit : took from delta_stepping_shortest_paths.hpp
  bool dist_changed = false;
  while (new_distance < old_dist) {
    last_old_dist = old_dist;
#ifdef CYCLE_TIME
    valstart = clock();
#else
    valstart = std::chrono::high_resolution_clock::now();
#endif
    old_dist = val_compare_and_swap(&distance_map[v], old_dist, new_distance);
#ifdef CYCLE_TIME
    valend = clock();
#else
    valend = std::chrono::high_resolution_clock::now();
#endif
    if (last_old_dist == old_dist) {
      // we got a better distance - insert it to work queue
#ifdef COLLECT_STATS
      if(old_dist < std::numeric_limits<Distance>::max()) {
	//	std::cout << "i:t=" << v << ":s=" << s << std::endl;
	++invalidated[amplusplus::detail::get_thread_id()];
      } else {
	++useful[amplusplus::detail::get_thread_id()];
      }
#endif
      dist_changed = true;
#ifdef CYCLE_TIME
      relstart = clock();
#else
      relstart = std::chrono::high_resolution_clock::now();
#endif
      BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
	Vertex u = target(e, graph);
#ifdef CYCLE_TIME
	ownstart = clock();
#else
	ownstart = std::chrono::high_resolution_clock::now();
#endif
	amplusplus::transport::rank_type dest = get(owner, u);
	last_owner = dest;
#ifdef CYCLE_TIME
	ownend = clock();
#else
	ownend = std::chrono::high_resolution_clock::now();	    
#endif

#ifdef NO_COALESCING
	work_item_t* pwi = new work_item_t(u, std::make_pair(std::make_pair(v, new_distance + 1), get_current_time()));
	mt.message_being_built(dest);
	//	msg_types[tid]->message_being_built(dest);
	//	msg_types[tid]->send(pwi, 1, dest, empty_deleter(pwi));
#ifdef CYCLE_TIME
	sendstart = clock();
#else
	sendstart = std::chrono::high_resolution_clock::now();
#endif
	mt.send(pwi, 1, dest, empty_deleter(pwi));
#ifdef CYCLE_TIME
	sendend = clock();
#else
	sendend = std::chrono::high_resolution_clock::now();
	totalsend += std::chrono::duration_cast<std::chrono::nanoseconds>(sendend - sendstart).count();
	++totalsendcnt;
#endif
#else
	work_item_t wi(u, std::make_pair(v, new_distance + 1));
	relax_msg.send(wi);
#endif
      }
#ifdef CYCLE_TIME
      relend = clock();
#else
      relend = std::chrono::high_resolution_clock::now();
#endif
 
      break;
    }
  }

#ifdef CYCLE_TIME
  auto end = clock();
#else
  auto end = std::chrono::high_resolution_clock::now();
#endif

  if (transport.rank() == 0) {
#ifdef CYCLE_TIME
    if (dist_changed) {
      std::cout << "Execution time : " 
		<< (end-start) << " (clicks), distance changed ? true "
		<< " Compared and swap time : " << (valend-valstart)
		<< " Relation time : " << (relend-relstart)
		<< " Out degree : " << boost::out_degree(v, graph) 
		<< " Owner check time : " << (ownend - ownstart)
		<< " Send time : " << (sendend - sendstart)
		<< " Destination : " << last_owner
		<< std::endl;
    } else {
      std::cout << "Execution time : " 
		<< end-start << " (nsec), distance changed ? false" << std::endl;
    }
#else
    if (dist_changed) {

      unsigned int sendtime = 0;
      if (totalsendcnt != 0)
	sendtime = (totalsend / totalsendcnt);

      /*      std::cout << "Execution time : " << 
	std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << " (nsec), distance changed ? true "
		<< " Compared and swap time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(valend-valstart).count() 
		<< " Relation time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(relend-relstart).count() 
		<< " Out degree : " << boost::out_degree(v, graph) 
		<< " Owner check time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(ownend - ownstart).count()
		<< " Send time : " << sendtime
		<< " Send count : " << totalsendcnt
		<< " Destination : " << last_owner
		<< std::endl;*/
    } else {
      //std::cout << "Execution time : " 
      //		<< std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << " (nsec), distance changed ? false" << std::endl;
    }
#endif
  }
#ifdef COLLECT_STATS
  if (!dist_changed) { 
    ++rejected[amplusplus::detail::get_thread_id()];
  }
#endif
}

template<BFS_CHAOTIC_PARMS>
struct BFS_CHAOTIC_TYPE::processing_function {
  processing_function() : self(NULL) {}
  processing_function(bfs_chaotic& self) : self(&self) {}

#ifdef NO_COALESCING
  void operator() (int source, work_item_t* pdata, size_t count) const {
    int tid = amplusplus::detail::get_thread_id();
    self->relax(*pdata, tid);
  }
#else
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
    self->relax(data, tid);
  }
#endif

protected:
  bfs_chaotic* self;
};

} // end namespace distributed
} // end namespace graph
} // end namespace boost

#endif
