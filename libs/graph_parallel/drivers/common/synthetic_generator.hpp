// Copyright (C) 2017 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine
#ifndef DRIVER_SYNTHETIC_GEN
#define DRIVER_SYNTHETIC_GEN

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <inttypes.h>
#include <cstdlib>
#include <math.h>

#define AMPLUSPLUS_PRINT_HIT_RATES
#define BFS_SV_CC_HACK
#define DISABLE_BFS_BITMAP
#ifdef __bgp__
#define BGP_REPORT_MEMORY 1 // if >1 all ranks will report memory usage
#endif
#define BARRIER_TRANS
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
#include <boost/graph/distributed/selector.hpp>
#include <boost/graph/recursive_rmat_generator.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/graph500_generator.hpp>
#include <boost/graph/permute_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/relax.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>
#include <boost/accumulators/accumulators.hpp>

#include <cds/init.h>       // for cds::Initialize and cds::Terminate
#include <cds/gc/hp.h>      // for cds::HP (Hazard Pointer) SMR

#include <parallel/algorithm>
#include <algorithm> // for std::min, std::max
#include <functional>
#include <boost/graph/util/utils.hpp>

// Params class -- try to keep this simple
class graph_gen_params {
public:
  graph_gen_params(): scale(20),
		      n((unsigned long long)(1) << scale),
		      gen_weights(true),
		      max_weight(255),
		      distribution_type(block),
		      block_size(10),
		      num_threads(16),
		      rmata(INITIATOR_A_NUMERATOR),
		      rmatbc(INITIATOR_BC_NUMERATOR),
		      preprocess(false)
		      
  {}

  int scale;
  unsigned long long n; // number of vertices
  
  // edges weights
  bool gen_weights;
  weight_type max_weight;

  distribution_t distribution_type;
  size_t block_size;
  size_t distribution_coalescing_size = 1 << 17;

  double edge_factor = 16;
  double edge_list_reserve_factor = 1.15;

  int num_threads;

  uint64_t seed64 = 12345;

  std::vector<routing_type> routing = {rt_none};

  // RMAT parameters -- These values will be divided by 10000
  int32_t rmata;
  int32_t rmatbc;

  // Preprocessing
  bool preprocess;

  void print() {
    assert(_RANK != -1);
    if (_RANK == 0) {
      std::cout << "============= Printing Graph Generation Parameters =============" << std::endl;
      std::cout << "Scale = " << scale << std::endl;
      std::cout << "n = " << n << std::endl;
      std::cout << "Generate Weights = " << gen_weights << std::endl;
      std::cout << "Max weight = " << max_weight << std::endl;
      std::cout << "Distribution type = " << distribution_type << std::endl;
      std::cout << "Block size = " << block_size << std::endl;
      std::cout << "Distribution Coalescing = " << distribution_coalescing_size << std::endl;
      std::cout << "Edge factor = " << edge_factor << std::endl;
      std::cout << "Edge list reserve factor = " << edge_list_reserve_factor << std::endl;
      std::cout << "Threads = " << num_threads << std::endl;
      std::cout << "Seed = " << seed64 << std::endl;
      assert(routing.size() >= 1);
      std::cout << "Routing = " << routing[0] << std::endl;
      std::cout << "RMAT A = " << rmata << std::endl; 
      std::cout << "RMAT B = C = " << rmatbc << std::endl; 
      std::cout << "Preprocess : " << preprocess << std::endl;
      std::cout << "================================================================" << std::endl;
    }
  }

  std::string get_graph_type() {
    std::string graph_type = "unknown";
    if ((rmata == 5700) &&
	(rmatbc == 1900))
      graph_type = "RMAT1";

    if ((rmata == 5000) &&
	(rmatbc == 1000))
      graph_type = "RMAT2";

    return graph_type;
  }

  std::string get_generator_specific_input() {
    std::stringstream ss;
    ss << ", RMAT-A :" << rmata
       << ", RMAT-BC :" << rmatbc
       << ", SCALE : " << scale
       << ", DEGREE : " << edge_factor;

    return ss.str();
  }


  bool load_graph_file() {
    return false;
  }

  bool parse(int argc, char* argv[]) {

    bool graph_type_specified = false;

    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--scale") {
	scale = boost::lexical_cast<size_t>( argv[i+1] );
	n = (unsigned long long)(1) << scale;
      }

      if (arg == "--no-weight") {
	gen_weights = false;
      }

      if (arg == "--max-weight") {
	max_weight = boost::lexical_cast<weight_type>( argv[i+1] );
      }

      if (arg == "--distribution-type") {
	std::string distt = boost::lexical_cast<std::string>( argv[i+1] );
	if (distt == "block")
	  distribution_type = block;
	else if (distt == "cyclic")
	  distribution_type = cyclic;
	else
	  std::cerr << "Invalid distribution type. Available types are block and cyclic" << std::endl;
      }

      if (arg == "--block-size") {
	block_size = boost::lexical_cast<size_t>( argv[i+1] );
      }

      if (arg == "--distribution-coalescing-size") {
	distribution_coalescing_size = boost::lexical_cast<size_t>( argv[i+1] );
      }

      if (arg == "--preprocess")
	preprocess = true;

      if (arg == "--edge-factor") {
	edge_factor = boost::lexical_cast<double>( argv[i+1] );
      }

      if (arg == "--edge-list-reserve-factor") {
	edge_list_reserve_factor = boost::lexical_cast<double>( argv[i+1] );
      }

      if (arg == "--num-threads") {
	num_threads = boost::lexical_cast<int>( argv[i+1] );
      }

      if (arg == "--seed") {
	seed64 = boost::lexical_cast<uint64_t>( argv[i+1] );
      }

      if (arg == "--routing") {
	std::string rt = boost::lexical_cast<std::string>( argv[i+1] );
	if (rt == "rt_none")
	  routing.push_back(rt_none);
	else if (rt == "rt_rook")
	  routing.push_back(rt_rook);
	else if (rt == "rt_hypercube")
	  routing.push_back(rt_hypercube);
	else 
	  std::cerr << "Invalid routing type. Available routing types -- rt_none, rt_rook and rt_hypercube"
		   << std::endl;
      }

      if (arg == "--rmat1") {

	if (graph_type_specified) {
	  std::cerr << "Graph type already specified. Cannot specify again." << std::endl;
	  return false;
	}

	rmata = 5700;
	rmatbc = 1900;
	graph_type_specified = true;
      }

      if (arg == "--rmat2") {
	if (graph_type_specified) {
	  std::cerr << "Graph type already specified. Cannot specify again." << std::endl;
	  return false;
	}

	rmata = 5000;
	rmatbc = 1000;
	graph_type_specified = true;
      }

      if (arg == "--rmata") {
	if (graph_type_specified) {
	  std::cerr << "Graph type already specified. Cannot specify again." << std::endl;
	  return false;
	}

	rmata = boost::lexical_cast<int>( argv[i+1] );
	assert(rmata <= INITIATOR_DENOMINATOR);
      }

      if (arg == "--rmatbc") {
	if (graph_type_specified) {
	  std::cerr << "Graph type already specified. Cannot specify again." << std::endl;
	  return false;
	}

	rmatbc = boost::lexical_cast<int>( argv[i+1] );
	assert(rmatbc <= INITIATOR_DENOMINATOR);
      }

      assert((rmata + rmatbc*2) < INITIATOR_DENOMINATOR);
    }

    return true;
  }
};

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


// For distribution
template<typename MessageType>
class threaded_distributor {
public:
  threaded_distributor(MessageType& mtype) : msg_type(mtype) {}
  void operator()() {
  }
  template <typename InIter, typename Flip, typename SelfLoopDetector>
  void operator()(int tid, InIter begin, InIter end, const Flip& flip,
		  const SelfLoopDetector sloopdet,
                  amplusplus::transport& trans) {
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      if (0 == tid) {
        int nthreads = trans.get_nthreads();
        t_bar.reset(new amplusplus::detail::barrier(nthreads));
      }

      // for thread barrier to initialize
      { amplusplus::scoped_epoch epoch(trans); }
      
      //      debug("starting concurrent distribute", tid);
      concurrent_distribute(begin, end, flip, sloopdet, trans, msg_type);

      t_bar->wait();
    }
  }

private:
  MessageType& msg_type;
  boost::shared_ptr<amplusplus::detail::barrier> t_bar;
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

// For distribution -- to generate an undirected graph with weights
struct flip_pair_w_weight {
  template <typename T>
  T operator()(const T& x) const {
    // flipping edges by keeping same weight
    return std::make_pair(x.second.first, std::make_pair(x.first, x.second.second));
  }
};

// For distribution -- to generate an undirected graph
struct flip_pair_wo_weight {
  template <typename T>
  T operator()(const T& x) const {
    // flipping edges by keeping same weight
    return std::make_pair(x.second, x.first);
  }
};

// Exclude self loops for weighted edges
struct a_self_looop_w_weight {
  template <typename T>
  bool operator()(const T& x) const {
    // flipping edges by keeping same weight
    return (x.first == x.second.first);
  }
};

// Exclude self-loops for unweighted edges
struct a_self_loops_wo_weight {
  template <typename T>
  bool operator()(const T& x) const {
    // flipping edges by keeping same weight
    return (x.first == x.second);
  }
};


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


class synthetic_generator {
  typedef boost::compressed_sparse_row_graph<boost::directedS, boost::no_property, WeightedEdge, 
					     boost::no_property, boost::distributedS<unsigned long long> > Digraph;
  typedef boost::graph_traits<Digraph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Digraph>::edges_size_type edges_size_type;
  typedef boost::graph_traits<Digraph>::vertices_size_type vertices_size_type;
  typedef boost::property_map<Digraph, boost::vertex_index_t>::type VertexIndexMap;
  typedef boost::property_map<Digraph, boost::vertex_owner_t>::const_type OwnerMap;

public:
  synthetic_generator(amplusplus::transport& transport):
    trans(transport), g(NULL) {}

 template<typename Graph500Iter, 
	  typename Routing, 
	  typename Distribution, 
	  typename EdgeSize, 
	  typename UInt, 
	  typename Edges, 
	  typename WeightIterator>
 void do_distribute(const graph_gen_params& params,
		    Distribution &distrib,
		    EdgeSize e_start, 
		    EdgeSize e_count, 
		    UInt a, 
		    UInt b, 
		    Edges &edges, 
		    WeightIterator weight_iter) {

    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, Routing> MessageGenerator;
    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(params.distribution_coalescing_size
#ifdef NO_COALESCING
									    , true)
#else
									    )
#endif
      , Routing(trans.rank(), trans.size()));

    time_type start = get_time();

    // Multithreaded distribution
    trans.set_nthreads(params.num_threads);

#ifdef PRINT_DEBUG
    std::cout << "Thread start : " << e_start << "Total edge count : " << e_count
	      << " for rank :" << trans.rank()
	      << std::endl;
#endif

    Edges* thread_edges = new Edges[params.num_threads];
    EdgeSize std_count_per_thread = (e_count + params.num_threads - 1) / params.num_threads;
    //   EdgeSize remainder = e_count % num_threads;

#ifdef PRINT_DEBUG
    std::cout << "std_count_per_thread : " << std_count_per_thread << std::endl;
#endif

    EdgeSize thread_start = e_start;
    EdgeSize count_per_thread = std_count_per_thread;

    // create message type
    typedef typename std::iterator_traits<Graph500Iter>::value_type value_type;
    typedef boost::detail::vector_write_handler<Edges> vec_write_handler;

    // not sure we need this here
    amplusplus::register_mpi_datatype<value_type>();
    
    typedef typename MessageGenerator::template call_result<value_type, vec_write_handler,
							    boost::detail::distrib_to_pair_owner_map_t<Distribution>, 
							    amplusplus::no_reduction_t>::type 
      write_msg_type;

    write_msg_type write_msg(msg_gen, trans, boost::detail::distrib_to_pair_owner_map<Distribution>(distrib), 
			     amplusplus::no_reduction);

    // set the handler (must be in main thread)
    write_msg.set_handler(vec_write_handler(thread_edges));
    threaded_distributor<write_msg_type> td(write_msg);

    boost::scoped_array<boost::thread> threads(new boost::thread[params.num_threads - 1]);
    for (int i = 0; i < params.num_threads - 1; ++i) {
#ifdef PRINT_DEBUG
      if (_RANK == 0) {
	std::cout << "Thread : " << i+1 << " start : " << thread_start
		  << " thread end : " << (thread_start + count_per_thread)
		  << std::endl;
      }
#endif
      boost::thread thr(boost::ref(td), i+1, Graph500Iter(params.scale, thread_start, a, b, weight_iter, params.rmata, params.rmatbc),
			Graph500Iter(params.scale, thread_start + count_per_thread, a, b, weight_iter, params.rmata, params.rmatbc),
			flip_pair_w_weight(), a_self_looop_w_weight(), trans);
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

    td(0, Graph500Iter(params.scale, thread_start, a, b, weight_iter),
       Graph500Iter(params.scale, thread_start + count_per_thread, a, b, weight_iter),
       flip_pair_w_weight(), a_self_looop_w_weight(), trans);

    // Wait till all threads finish
    for (int i = 0; i < params.num_threads - 1; ++i)
      threads[i].join();

    // Merge all edges - in threads
    //    time_type m_start = get_time();

    typename Edges::size_type total_sz = 0;
    std::vector<typename Edges::size_type> positions(params.num_threads);
    for (int i = 0; i < params.num_threads; ++i) {
      positions[i] = total_sz;
      total_sz += thread_edges[i].size();
    }

    edges.resize(total_sz);
    threaded_merger<Edges> tm(edges);

    boost::scoped_array<boost::thread> m_threads(new boost::thread[params.num_threads - 1]);
    for (int i = 0; i < params.num_threads - 1; ++i) {
      typename Edges::iterator pos_ite = edges.begin();
      std::advance(pos_ite, positions[i+1]);
      boost::thread thr(boost::ref(tm), pos_ite, 
			thread_edges[i+1].begin(), thread_edges[i+1].end());
      m_threads[i].swap(thr);
    }

    tm(edges.begin(), thread_edges[0].begin(), thread_edges[0].end());

    // Wait till all threads finish
    for (int i = 0; i < params.num_threads - 1; ++i)
      m_threads[i].join();

    //time_type m_end = get_time();

    trans.set_nthreads(1);
    delete[] thread_edges;
    { amplusplus::scoped_epoch epoch(trans); }
    time_type end = get_time();
    if (_RANK==0)
      std::cout << "[INFO] Graph distribution permutation took " << (end-start)
		<< " time" << std::endl;

    // sort edges in parallel
    //__gnu_parallel::sort(edges.begin(), edges.end());
  }


  // distribution without weights -- need to find a better way to do this without duplicating code
 template<typename Graph500Iter, 
	  typename Routing, 
	  typename Distribution, 
	  typename EdgeSize, 
	  typename UInt, 
	  typename Edges>
 void do_distribute(const graph_gen_params& params,
		    Distribution &distrib,
		    EdgeSize e_start, 
		    EdgeSize e_count, 
		    UInt a, 
		    UInt b, 
		    Edges &edges) {

    typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, Routing> MessageGenerator;
    MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(params.distribution_coalescing_size
#ifdef NO_COALESCING
									    , true)
#else
									    )
#endif
      , Routing(trans.rank(), trans.size()));


    time_type start = get_time();

    // Multithreaded distribution
    trans.set_nthreads(params.num_threads);

#ifdef PRINT_DEBUG
    std::cout << "Thread start : " << e_start << "Total edge count : " << e_count
	      << " for rank :" << trans.rank()
	      << std::endl;
#endif

    Edges* thread_edges = new Edges[params.num_threads];
    EdgeSize std_count_per_thread = (e_count + params.num_threads - 1) / params.num_threads;
    //   EdgeSize remainder = e_count % num_threads;

#ifdef PRINT_DEBUG
    std::cout << "std_count_per_thread : " << std_count_per_thread << std::endl;
#endif

    EdgeSize thread_start = e_start;
    EdgeSize count_per_thread = std_count_per_thread;

    // create message type
    typedef typename std::iterator_traits<Graph500Iter>::value_type value_type;
    typedef boost::detail::vector_write_handler<Edges> vec_write_handler;

    // not sure we need this here
    amplusplus::register_mpi_datatype<value_type>();
    
    typedef typename MessageGenerator::template call_result<value_type, vec_write_handler,
							    boost::detail::distrib_to_pair_owner_map_t<Distribution>, 
							    amplusplus::no_reduction_t>::type 
      write_msg_type;

    write_msg_type write_msg(msg_gen, trans, boost::detail::distrib_to_pair_owner_map<Distribution>(distrib), 
			     amplusplus::no_reduction);

    // set the handler (must be in main thread)
    write_msg.set_handler(vec_write_handler(thread_edges));
    threaded_distributor<write_msg_type> td(write_msg);

    boost::scoped_array<boost::thread> threads(new boost::thread[params.num_threads - 1]);
    for (int i = 0; i < params.num_threads - 1; ++i) {
#ifdef PRINT_DEBUG
      if (_RANK == 0) {
	std::cout << "Thread : " << i+1 << " start : " << thread_start
		  << " thread end : " << (thread_start + count_per_thread)
		  << std::endl;
      }
#endif
      boost::thread thr(boost::ref(td), i+1, Graph500Iter(params.scale, thread_start, a, b, params.rmata, params.rmatbc),
			Graph500Iter(params.scale, thread_start + count_per_thread, a, b, params.rmata, params.rmatbc),
			flip_pair_wo_weight(), a_self_loops_wo_weight(), trans);
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

    td(0, Graph500Iter(params.scale, thread_start, a, b),
       Graph500Iter(params.scale, thread_start + count_per_thread, a, b),
       flip_pair_wo_weight(), a_self_loops_wo_weight(), trans);

    // Wait till all threads finish
    for (int i = 0; i < params.num_threads - 1; ++i)
      threads[i].join();

    // Merge all edges - in threads
    //    time_type m_start = get_time();

    typename Edges::size_type total_sz = 0;
    std::vector<typename Edges::size_type> positions(params.num_threads);
    for (int i = 0; i < params.num_threads; ++i) {
      positions[i] = total_sz;
      total_sz += thread_edges[i].size();
    }

    edges.resize(total_sz);
    threaded_merger<Edges> tm(edges);

    boost::scoped_array<boost::thread> m_threads(new boost::thread[params.num_threads - 1]);
    for (int i = 0; i < params.num_threads - 1; ++i) {
      typename Edges::iterator pos_ite = edges.begin();
      std::advance(pos_ite, positions[i+1]);
      boost::thread thr(boost::ref(tm), pos_ite, 
			thread_edges[i+1].begin(), thread_edges[i+1].end());
      m_threads[i].swap(thr);
    }

    tm(edges.begin(), thread_edges[0].begin(), thread_edges[0].end());

    // Wait till all threads finish
    for (int i = 0; i < params.num_threads - 1; ++i)
      m_threads[i].join();

    //time_type m_end = get_time();

    trans.set_nthreads(1);
    delete[] thread_edges;
    { amplusplus::scoped_epoch epoch(trans); }
    time_type end = get_time();
    if (_RANK == 0)
      std::cout << "[INFO] Graph distribution permutation took " << (end-start)
		<< " time" << std::endl;

    // sort edges in parallel
    //__gnu_parallel::sort(edges.begin(), edges.end());
  }



  struct vertex_comparator {
    template<typename Edge>
    bool operator()(const Edge& e1, const Edge& e2) {
      if (e1.first < e2.first)
	return true;
      else if (e1.first == e2.first) {
	return e1.second < e2.second;
      } else
	return false;
    }
  };

  // Preprocess the graph
  // 1. Remove self loops
  // 2. Remove parallel edges
  template<typename edge_t>
  void preprocess(std::vector<edge_t>& edges) {
    std::sort(edges.begin(), edges.end(), vertex_comparator());
    edges.erase(std::unique(edges.begin(), edges.end(), 
			    [](edge_t e1, edge_t e2) { 
			      return ((e1.first == e2.first) && (e1.second == e2.second));
			    }), edges.end());

#ifdef PRINT_DEBUG
    for (auto i=0; i < edges.size(); ++i) {
      std::cout << "(" << edges[i].first << ", " << edges[i].second << ")" << std::endl;
    }
#endif
  }


  void generate_graph_w_weights(const graph_gen_params& params) {

    if (_RANK == 0)
      std::cout << "[INFO] Generating graph with edges" << std::endl;

    // Seed general-purpose RNG
    boost::rand48 gen, synch_gen;
    gen.seed(params.seed64);
    synch_gen.seed(params.seed64);

    typedef generator_iterator<boost::uniform_int<weight_type>, boost::rand48> weight_iterator_t;
    weight_iterator_t gi 
      = make_generator_iterator(gen, boost::uniform_int<weight_type>(1, params.max_weight));


    typedef std::vector<std::pair<edges_size_type, std::pair<edges_size_type, weight_type> > > edge_with_weight_t;
    edge_with_weight_t edges;
    boost::parallel::variant_distribution<vertices_size_type> distrib;

    edges_size_type m = static_cast<edges_size_type>(floor(params.n * params.edge_factor));
    edges.reserve(static_cast<edges_size_type>(floor(params.edge_list_reserve_factor * 2 * m / trans.size())));

    if (params.distribution_type == cyclic) {
      distrib = boost::parallel::oned_block_cyclic<vertices_size_type>(trans, params.block_size);
    } else //(distribution_type == block)
      distrib = boost::parallel::block<vertices_size_type>(trans, params.n);


    boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());

    edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
    edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

    // Permute and redistribute copy constructs the input iterator
    uint64_t a = rand_64(gen);
    uint64_t b = rand_64(gen);

    // Build a graph to test with
    typedef boost::graph500_iterator<Digraph, 
				     generator_iterator<boost::uniform_int<weight_type>, 
							boost::rand48>, weight_type > Graph500Iter;
    assert(!params.routing.empty());

    // As for now, we assume that only the first routing from the list is used to generate the graph.  
    // It seems that it does not really matter which routing we use here since generation of the graph is not a part of the performance test.
    if (params.routing[0] == rt_none) {
      do_distribute<Graph500Iter, amplusplus::no_routing>(params, distrib, e_start, e_count, a, b, edges, gi);
    } else if (params.routing[0] == rt_hypercube) {
      do_distribute<Graph500Iter, amplusplus::hypercube_routing>(params, distrib, e_start, e_count, a, b, edges, gi);
    } else if (params.routing[0] == rt_rook) {
      do_distribute<Graph500Iter, amplusplus::rook_routing>(params, distrib, e_start, e_count, a, b, edges, gi);
    }

    typedef select_edge<edges_size_type, weight_type> EdgeSelectFunction;
    typedef boost::transform_iterator<EdgeSelectFunction, edge_with_weight_t::iterator> edge_only_iterator;

    typedef select_edge_weight<edges_size_type, weight_type> WeightSelectFunction;
    typedef boost::transform_iterator<WeightSelectFunction, edge_with_weight_t::iterator> weight_only_iterator;

    edge_only_iterator edge_begin(edges.begin(), EdgeSelectFunction()), 
      edge_end(edges.end(), EdgeSelectFunction());

    weight_only_iterator weight_begin(edges.begin(), WeightSelectFunction());

    time_type gt1 = get_time();
    g = new Digraph(boost::edges_are_unsorted_multi_pass, edge_begin, edge_end,
		    weight_begin, params.n, trans, distrib);
    time_type gt2 = get_time();

    if (_RANK == 0)
      std::cout << "[INFO] Local graph creation time : " << (gt2-gt1)
		<< std::endl;
    // Clear edge array above
    edges.clear();

  }

  void generate_graph_wo_weights(const graph_gen_params& params) {

    if (_RANK == 0)
      std::cout << "[INFO] Generating graph with out edges" << std::endl;

    // Seed general-purpose RNG
    boost::rand48 gen, synch_gen;
    gen.seed(params.seed64);
    synch_gen.seed(params.seed64);

    typedef std::pair<edges_size_type, edges_size_type> edge_t;
    typedef std::vector<edge_t> edges_t;
    edges_t edges;
    boost::parallel::variant_distribution<vertices_size_type> distrib;

    edges_size_type m = static_cast<edges_size_type>(floor(params.n * params.edge_factor));
    edges.reserve(static_cast<edges_size_type>(floor(params.edge_list_reserve_factor * 2 * m / trans.size())));

    if (params.distribution_type == cyclic) {
      distrib = boost::parallel::oned_block_cyclic<vertices_size_type>(trans, params.block_size);
    } else //(distribution_type == block)
      distrib = boost::parallel::block<vertices_size_type>(trans, params.n);


    boost::uniform_int<uint64_t> rand_64(0, std::numeric_limits<uint64_t>::max());

    edges_size_type e_start = trans.rank() * (m + trans.size() - 1) / trans.size();
    edges_size_type e_count = (std::min)((m + trans.size() - 1) / trans.size(), m - e_start);

    // Permute and redistribute copy constructs the input iterator
    uint64_t a = rand_64(gen);
    uint64_t b = rand_64(gen);

    // Build a graph to test with
    typedef boost::graph500_iterator<Digraph> Graph500Iter;
    assert(!params.routing.empty());

    // As for now, we assume that only the first routing from the list is used to generate the graph.  
    // It seems that it does not really matter which routing we use here since generation of the graph is not a part of the performance test.
    if (params.routing[0] == rt_none) {
      do_distribute<Graph500Iter, amplusplus::no_routing>(params, distrib, e_start, e_count, a, b, edges);
    } else if (params.routing[0] == rt_hypercube) {
      do_distribute<Graph500Iter, amplusplus::hypercube_routing>(params, distrib, e_start, e_count, a, b, edges);
    } else if (params.routing[0] == rt_rook) {
      do_distribute<Graph500Iter, amplusplus::rook_routing>(params, distrib, e_start, e_count, a, b, edges);
    }

    if (params.preprocess) {
      if (_RANK == 0)
	std::cout << "[INFO] Preprocessing the graph -- removing parallel edges" << std::endl;

      time_type pt1 = get_time();
      preprocess<edge_t>(edges);
      time_type pt2 = get_time();

      if (_RANK == 0)
	std::cout << "[INFO] The preprocessing time : " << pt2-pt1 << std::endl;
    }


    time_type gt1 = get_time();
    g = new Digraph(boost::edges_are_unsorted_multi_pass, edges.begin(), edges.end(),
	      params.n, trans, distrib);
    time_type gt2 = get_time();

    if (_RANK == 0)
      std::cout << "[INFO] Local graph creation time : " << (gt2-gt1)
	      << std::endl;
    // Clear edge array above
    edges.clear();
  }

  void generate_graph(const graph_gen_params& params) {
    time_type start = get_time();

    if (params.gen_weights)
      generate_graph_w_weights(params);
    else
      generate_graph_wo_weights(params);

    time_type end = get_time();
    if (_RANK == 0) 
      std::cout << "[INFO] Graph generation took " << end - start 
	      << "(sec) time." << std::endl;
  }

  Digraph* graph() {
    return g;
  }

  // MUST call finalize release memory associated with the graph
  void finalize() {
    delete g;
  }

private:
  amplusplus::transport& trans;
  Digraph* g;

};

#endif

