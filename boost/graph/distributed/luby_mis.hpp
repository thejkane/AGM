// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Maximal Independent Set Algortihms================//
// This implements Luby algorithms: 1. A - Luby algorithm A 
// MIS algorithm, 2. AV1 - A variation of Luby's algorithm A,
// 3. AV2 - Another variation of Luby A, 4. Luby B algorithm.  
// ordering (use MIS_PRIORITY preprocessor macro)
//===========================================================//

#ifndef BOOST_GRAPH_LUBY_MIS
#define BOOST_GRAPH_LUBY_MIS

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
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/random_device.hpp>

//#include <pat_api.h>

#define REMOVE_VERTEX 1
#define REMOVE_REPLY 2
#define UNFIX_VERTEX 3
#define IN_MIS_QUERY 4
#define IN_MIS_REPLY 5


typedef int action_t;
typedef int mis_thread_t;

#if __x86_64__ || __ppc64__
typedef int64_t random_t;
#else
#warning("Compiling in 32 bits")
typedef int32_t random_t;
#endif

namespace boost { namespace graph { namespace distributed {

template<typename SelectFunctor>
struct processing_function;

struct minimum_pair_first {
  template<typename T>
  const T& operator()(const T& x, const T& y) const { return x.first < y.first ? x : y; }

  template<typename F>
  struct result {
    typedef typename boost::function_traits<F>::arg1_type type;
  };
};


//==========================================================//
//    Threaded Set
//==========================================================//
template<typename T>
class multi_thread_set {

public:
  multi_thread_set(int threads) : nthreads(threads) {
    threaded_set.resize(nthreads);
  }

  bool insert(T v, int tid) {
    assert(tid < nthreads);
    threaded_set[tid].insert(v);
  }

  size_t size(int tid) {
    return threaded_set[tid].size();
  }

  size_t size() {
    size_t sz = 0;
    for (int i=0; i < nthreads; ++i) {
      sz += threaded_set[i].size();
    }
    return sz;
  }

  void clear() {
    for(int i=0; i < nthreads; ++i) {
      threaded_set[i].clear();
    }
  }

  bool exists(T v, int tid) {
    // first check in current thread
    if (threaded_set[tid].find(v) != threaded_set[tid].end()) {
      // found in the current thread set
      return true;
    } else {
      // go through other set and see whether vertex exists
      for (int i=0; i < nthreads; ++i) {
	if (i != tid) {
	  if (threaded_set[i].find(v) != threaded_set[i].end()) { 
	    return true;
	  }
	}
      }
    }

    // did not find the vertex in any of the sets
    return false;
  }

private:
  int nthreads;
  std::vector<std::set<T> > threaded_set;

};


//======================================Select A Functor ======================================//
// Functor for Select A
template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class SelectA {

  // Following definitions must match with the definitions in the graph
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef append_buffer<Vertex> Bucket_t;
  typedef std::pair<Vertex, std::pair<Vertex, std::pair<random_t, action_t> > > work_item_t;
  //  typedef std::set<Vertex> DeletedVertices;
  typedef multi_thread_set<Vertex> DeletedVertices;
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef SelectA<Graph, MISMap, RandomMap, MessageGenerator> self_type;

  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function<self_type>, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;


public:
  SelectA(Graph& graph, 
	  MISMap& maxmap,
	  RandomMap& random,
	  amplusplus::transport& t) : g(graph), p_local_mis(NULL), 
				      p_owner(NULL),
				      mis(maxmap), rmap(random),
				      p_deleted_vertices(NULL),
				      p_relax_msg(NULL),
				      transport(t), first_call(true){}

  void initialize(Bucket_t* lmis,
		  OwnerMap* owner,
		  DeletedVertices* delv,
		  RelaxMessage* relmsg,
		  unsigned long long numv,
		  int nthreads) {

    p_local_mis = lmis;
    p_owner = owner;
    p_deleted_vertices = delv;
    p_relax_msg = relmsg;

#ifdef MIS_STATS
    //    selecttime.resize(nthreads);
    //mergetime.resize(nthreads);
    rndcalctime.resize(nthreads);
    subgcomptime.resize(nthreads);
    iscalctime.resize(nthreads);
    localsetsz.resize(nthreads);
    removesetsz.resize(nthreads);
#endif

  }


#ifdef MIS_STATS
  void print_stats(int nthreads) {
    auto iterations = rndcalctime[0].size();
    std::cout << "Total number of iterations : " << iterations << std::endl;
    for (int i=0; i < nthreads; ++i) {
      for (int k=0; k < iterations; ++k) {
	//std::cout << "tid : " << i << ", ite : " << k << ", Random value calculate Time : " << rndcalctime[i][k] << std::endl;
	//std::cout << "tid : " << i << ", ite : " << k << ", Sub-graph compute Time : " << subgcomptime[i][k] << std::endl;
	std::cout << "tid : " << i << ", ite : " << k << ", Independent Set Calculation Time : " << iscalctime[i][k] << std::endl;
	std::cout << "tid : " << i << ", ite : " << k << ", Local Independent Set Size : " << localsetsz[i][k] << std::endl;
	std::cout << "tid : " << i << ", ite : " << k << ", Remove Set Size : " << removesetsz[i][k] << std::endl;
      }
    }

    std::cout << "========================== Aggregated Independent Set Sizes ==============================" << std::endl;
    //aggregate over threads
    std::vector<unsigned long> aggrlocalset;
    std::vector<unsigned long> aggrremoveset;
    aggrlocalset.resize(iterations, 0);
    aggrremoveset.resize(iterations, 0);

    for (int k=0; k < iterations; ++k) {
      for (int i=0; i < nthreads; ++i) {
	aggrlocalset[k] += localsetsz[i][k];
	aggrremoveset[k] += removesetsz[i][k];
      } 
    }

    for (int i=0; i < iterations; ++i) {
      std::cout << "ite : " << i << ", Assumed Independent Set Size : " << aggrlocalset[i] << std::endl;
      std::cout << "ite : " << i << ", Deleted Set Size : " << aggrremoveset[i] << std::endl;
    }

  }
#endif

  random_t power4(std::size_t numv) {
    int power = 4;
    random_t pn = std::pow(numv, power);

    // adjusting pn as necessary
    while ((pn <= 0) && (power > 0)) {
      --power;
      pn = std::pow(numv, power);
    }

#ifdef PRINT_DEBUG
    std::cout << "numv :" << numv << "Power value : " << power << ", pn : " << pn << std::endl;
#endif

    assert((pn > 0) && "n^4 exceeds maximum limit for unsigned long long");
    return pn;
  }

  void operator() (int tid, int nthreads,
		   std::size_t& unfix_set_sz,
		   shared_ptr<amplusplus::detail::barrier>& t_bar) {

    assert(p_local_mis != NULL);
    assert(p_deleted_vertices != NULL);
    assert(p_owner != NULL);
    assert(p_relax_msg != NULL);

    if (tid == 0) {
      if (first_call) {
	generators.resize(nthreads);

	for(int k=0; k < nthreads; ++k) {
	  auto seedv = k + (nthreads * transport.rank());
	  generators[k].seed(seedv);
	}

	first_call = false;
      }
    }
    // wait till above code is executed by main thread in every rank
    t_bar->wait();

#ifdef MIS_STATS
    time_type asubst = get_time();
#endif

    std::size_t all_sizes = 0;
    {
      amplusplus::scoped_epoch_value epoch(transport, unfix_set_sz, all_sizes);
      BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
	if (mis[u] == MIS_UNFIX) {
	  ++unfix_set_sz;
	}
      }
    }

#ifdef MIS_STATS
    time_type esubst = get_time();
    subgcomptime[tid].push_back(esubst - asubst);
#endif

#ifdef MIS_STATS
  time_type srndtm = get_time();
#endif
    random_t pn = power4(all_sizes);
    boost::random::uniform_int_distribution<random_t> randq(1, pn);

    BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
      if (mis[u] == MIS_UNFIX) {
	rmap[u] = randq(generators[tid]);
	add_to_local_set(u);
      }
    }

#ifdef MIS_STATS
  time_type erndtm = get_time();
  rndcalctime[tid].push_back(erndtm-srndtm);
  localsetsz[tid].push_back(p_local_mis->size());
#endif


    { amplusplus::scoped_epoch epoch(transport); }

#ifdef MIS_STATS
    time_type t0 = get_time();
#endif

    { 
#ifdef PRINT_DEBUG
      if ((tid == 0) && (transport.rank() == 0))
	std::cout << "Local MIS size : " << p_local_mis->size() << std::endl;
#endif

      amplusplus::scoped_epoch epoch(transport); 

      for (typename Bucket_t::size_type i = tid ; i < p_local_mis->size() ; i += nthreads) {
	Vertex v = (*p_local_mis)[i];
	// synchronize all the nodes
	// thread barrier is impicit ?
	//    { amplusplus::scoped_epoch epoch(transport); }

	// There can be many edges between two vertices
	// For algorithm to work we need to make sure
	// each adjacen vertex is counted only once
	// we will have adjacencies set per each thread
	std::set<Vertex> adjacencies;

	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  if (v != u) { // ignore self-loops
	    if (adjacencies.insert(u).second) {
	
	      // is u owned by current rank ? then fine, proceed with the algorithm
	      // otherwise send a message to owner with continuation data
	      if (transport.rank() == get(*p_owner, u)) {
		if (mis[u] == MIS_UNFIX) { // This is necessary. Because we are working on the subgraph.
#ifdef PRINT_DEBUG
		  std::cout << "r-v : " << rmap[v] << ", r-u : " << rmap[u] << " local : " 
			    << (transport.rank() == get(*p_owner, u)) << std::endl;
#endif
		  if (rmap[v] > rmap[u]) {
		    // souce must be local to the rank
		    remove_from_local_set(v, tid);
#ifdef PRINT_DEBUG
		    std::cout << "Removing Vertex-a " << v << "- (" << v << "," << u << ")" << std::endl;
#endif
		  } else if (rmap[v] == rmap[u]){
		    // TODO:use vertex id to break the symmetry -- Not exactly Luby. Need to discuss with Marcin
		    if (v < u) {
		      remove_from_local_set(v, tid);
#ifdef PRINT_DEBUG
		      std::cout << "Removing Vertex-c " << v << "- (" << v << "," << u << ")" << std::endl;
#endif
		    }
	  
		  } else {
		    remove_from_local_set(u, tid);
#ifdef PRINT_DEBUG
		    std::cout << "Removing Vertex-b " << u << "- (" << v << "," << u << ")" << std::endl;
#endif
		  }
		}
	      } else {
		// u is in a remote rank
		// we need to send v, rmap[v], u, action
		// if rmap[v] > rmap[u], the remote rank needs to send a reply back to current rank : handler
		// if rmap[v] == rmap[u] and v < u remote rank needs to send a reply back to current rank : handler
		// if rmap[v] < rmap[u], the remote rank needs to remove u from its local set, no need to send a reply
		// work_item_t = < target, source, source_random, action >
		work_item_t wi = construct_wi(u, v, rmap[v], REMOVE_VERTEX);
		p_relax_msg->send(wi);
	      }
	    }
	  }
	}
      }
    }

#ifdef MIS_STATS
    time_type t1 = get_time();
    iscalctime[tid].push_back(t1-t0);
    removesetsz[tid].push_back(p_deleted_vertices->size());
#endif

  } // end of operator


  // Active Messages are handled in this function
  void process(const work_item_t& data, int tid) {
    // TODO : simplify message structure after we confirm
    // everything is working. I am suspecting we may need to
    // send a reply after removing vertex from local set
    Vertex u = data.first;
    Vertex v = data.second.first;
    random_t rmap_v = data.second.second.first;
    action_t action = data.second.second.second;

    if (action == REMOVE_VERTEX) {
      if (mis[u] == MIS_UNFIX) {
	if (rmap_v > rmap[u]) {
	  // send a reply -- TODO : We may can simplify the reply
	  work_item_t wi = construct_wi(v, u, rmap[u], REMOVE_REPLY);
	  p_relax_msg->send(wi);
	} else if (rmap_v == rmap[u]) {
	  if (v < u) {
	    // send a reply
	    work_item_t wi = construct_wi(v, u, rmap[u], REMOVE_REPLY);
	    p_relax_msg->send(wi);
	  } else
	    remove_from_local_set(u, tid);
	} else {
	  remove_from_local_set(u, tid);
	}
      }
    } else if (action == REMOVE_REPLY) {
      // just remove the vertex from local set
      remove_from_local_set(v, tid);
    } else if (action == UNFIX_VERTEX) {
      int expected = MIS_UNFIX;
      if (__atomic_compare_exchange_n(&mis[u], &expected, MIS_FIX0, false /*no weak*/, 
				      __ATOMIC_SEQ_CST, __ATOMIC_RELAXED)) {
      }

    } else
      assert(false);
  }

private:
  void add_to_local_set(Vertex v) {
    p_local_mis->push_back(v);
  }

  void remove_from_local_set(Vertex v, int tid) {
    // TODO this could be a map. Not sure which one is going to give
    // better performance. Experiment.
    //boost::mutex::scoped_lock
    //  lock(set_mutex);
    p_deleted_vertices->insert(v, tid);
  }

  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   random_t source_random,
			   action_t action) {
    work_item_t wi(d, std::make_pair(s, std::make_pair(source_random, action)));
    return wi;
  }


  //  random_t pn; // number of total vertices
  amplusplus::transport& transport;
  Graph& g;
  Bucket_t* p_local_mis;
  OwnerMap* p_owner;
  MISMap& mis;
  RandomMap& rmap;
  DeletedVertices* p_deleted_vertices;
  RelaxMessage* p_relax_msg;
  bool first_call;
  std::vector<boost::random::mt19937> generators;
  boost::mutex set_mutex;
  boost::mutex gen_mutex;

#ifdef MIS_STATS
  std::vector<std::vector<time_t> > rndcalctime;
  std::vector<std::vector<time_t> > mergetime;
  std::vector<std::vector<time_t> > subgcomptime;
  std::vector<std::vector<time_t> > iscalctime;
  std::vector<std::vector<unsigned long> > localsetsz;
  std::vector<std::vector<unsigned long> > removesetsz;
#endif


};

// SelectA functor generator.
struct select_a_functor_gen {
  template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
  struct select {
    typedef SelectA<Graph, MISMap, RandomMap, MessageGenerator> type;
  };
};

//====================== End of SelectA Functor ==================================//


      //%%%%%%%%%%%%%%%%%%%%%%%%%//

//======================================Select A-V2 Functor ======================================//
// Functor for Select AV2 
template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class SelectAV2 {

  // Following definitions must match with the definitions in the graph
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef append_buffer<Vertex> Bucket_t;
  typedef std::pair<Vertex, std::pair<Vertex, std::pair<random_t, action_t> > > work_item_t;
  //  typedef std::set<Vertex> DeletedVertices;
  typedef multi_thread_set<Vertex> DeletedVertices;
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef SelectAV2<Graph, MISMap, RandomMap, MessageGenerator> self_type;

  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function<self_type>, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;


public:
  SelectAV2(Graph& graph, 
	  MISMap& maxmap,
	  RandomMap& random,
	  amplusplus::transport& t) : g(graph), p_local_mis(NULL), 
				      p_owner(NULL),
				      mis(maxmap), rmap(random),
				      p_deleted_vertices(NULL),
				      p_relax_msg(NULL),
				      pn(0), transport(t), gen(std::time(0) + std::time(0) * t.rank()){}

  void initialize(Bucket_t* lmis,
		  OwnerMap* owner,
		  DeletedVertices* delv,
		  RelaxMessage* relmsg,
		  unsigned long long numv,
		  int nthreads) {

    p_local_mis = lmis;
    p_owner = owner;
    p_deleted_vertices = delv;
    p_relax_msg = relmsg;


    // checking the overflow of numbers
    // NOTE -- this only works with GCC (maybe clang too)
    random_t result1;
    random_t result2;
    if (!__builtin_mul_overflow(numv, numv, &result1)) {
      // not overflowing
      if (!__builtin_mul_overflow(result1, numv, &result2)) {
	// still not overflowing
	if (!__builtin_mul_overflow(result2, numv, &pn)) {
	  // not overflowing
	  std::cout << "successfully calculated n^4 :" << pn << std::endl;
	} else {
	  pn = std::numeric_limits<random_t>::max();
	  std::cout << "n^4 overflowing. settling for random_t MAX : " << pn << std::endl;
	}
      } else {
	pn = std::numeric_limits<random_t>::max();
	std::cout << "n^3 overflowing. settling for random_t  MAX : " << pn << std::endl;
      }
    } else {
      pn = std::numeric_limits<random_t>::max();
      std::cout << "n^2 overflowing. settling for random_t MAX : " << pn << std::endl;
    }


    //#ifdef PRINT_DEBUG
    std::cout << "numv :" << numv << ", pn : " << pn << std::endl;
    //#endif

    assert((pn > 0) && "n^4 exceeds maximum limit for unsigned long long");
  }

#ifdef MIS_STATS
  void print_stats(int nthreads) {
    /*    auto iterations = rndcalctime[0].size();
    std::cout << "Total number of iterations : " << iterations << std::endl;
    for (int i=0; i < nthreads; ++i) {
      for (int k=0; k < iterations; ++k) {
	std::cout << "Random value calculate Time : " << rndcalctime[i][k] << std::endl;
	std::cout << "Sub-graph compute Time : " << subgcomptime[i][k] << std::endl;
	std::cout << "Independent Set Calculation Time : " << iscalctime[i][k] << std::endl;
	std::cout << "Independent Set Size : " << issize[i][k] << std::endl;
      }
      }*/
  }
#endif


  void operator() (int tid, int nthreads,
		   std::size_t& unfix_set_sz,
		   shared_ptr<amplusplus::detail::barrier>& t_bar) {

    assert(p_local_mis != NULL);
    assert(p_deleted_vertices != NULL);
    assert(p_owner != NULL);
    assert(p_relax_msg != NULL);
    assert(pn != 0);

#ifdef PRINT_DEBUG
    std::cout << "The pn value : " << pn << std::endl;
#endif
    boost::random::uniform_int_distribution<random_t> randq(1, pn);

    BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
      if (get(*p_owner, u) == transport.rank()) {
	if (mis[u] == MIS_UNFIX) {
	  ++unfix_set_sz;
	  {
	    boost::mutex::scoped_lock
	      lock(gen_mutex);
	    rmap[u] = randq(gen);
	  }

	  add_to_local_set(u);
	}
      }
    }

#ifdef PRINT_DEBUG
    time_type t0 = get_time();
#endif

    { 
#ifdef PRINT_DEBUG
      if ((tid == 0) && (transport.rank() == 0))
	std::cout << "Local MIS size : " << p_local_mis->size() << std::endl;
#endif

      amplusplus::scoped_epoch epoch(transport); 

      for (typename Bucket_t::size_type i = tid ; i < p_local_mis->size() ; i += nthreads) {
	Vertex v = (*p_local_mis)[i];
	// synchronize all the nodes
	// thread barrier is impicit ?
	//    { amplusplus::scoped_epoch epoch(transport); }

	// There can be many edges between two vertices
	// For algorithm to work we need to make sure
	// each adjacen vertex is counted only once
	// we will have adjacencies set per each thread
	std::set<Vertex> adjacencies;

	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  if (v != u) { // ignore self-loops
	    if (adjacencies.insert(u).second) {
	
	      // is u owned by current rank ? then fine, proceed with the algorithm
	      // otherwise send a message to owner with continuation data
	      if (transport.rank() == get(*p_owner, u)) {
		if (mis[u] == MIS_UNFIX) { // This is necessary. Because we are working on the subgraph.
#ifdef PRINT_DEBUG
		  std::cout << "r-v : " << rmap[v] << ", r-u : " << rmap[u] << " local : " 
			    << (transport.rank() == get(*p_owner, u)) << std::endl;
#endif
		  if (rmap[v] > rmap[u]) {
		    // souce must be local to the rank
		    remove_from_local_set(v, tid);
#ifdef PRINT_DEBUG
		    std::cout << "Removing Vertex-a " << v << "- (" << v << "," << u << ")" << std::endl;
#endif
		  } else if (rmap[v] == rmap[u]){
		    // TODO:use vertex id to break the symmetry -- Not exactly Luby. Need to discuss with Marcin
		    if (v < u) {
		      remove_from_local_set(v, tid);
#ifdef PRINT_DEBUG
		      std::cout << "Removing Vertex-c " << v << "- (" << v << "," << u << ")" << std::endl;
#endif
		    }
	  
		  } else {
		    remove_from_local_set(u, tid);
#ifdef PRINT_DEBUG
		    std::cout << "Removing Vertex-b " << u << "- (" << v << "," << u << ")" << std::endl;
#endif
		  }
		}
	      } else {
		// u is in a remote rank
		// we need to send v, rmap[v], u, action
		// if rmap[v] > rmap[u], the remote rank needs to send a reply back to current rank : handler
		// if rmap[v] == rmap[u] and v < u remote rank needs to send a reply back to current rank : handler
		// if rmap[v] < rmap[u], the remote rank needs to remove u from its local set, no need to send a reply
		// work_item_t = < target, source, source_random, action >
		work_item_t wi = construct_wi(u, v, rmap[v], REMOVE_VERTEX);
		p_relax_msg->send(wi);
	      }
	    }
	  }
	}
      }
    }

#ifdef PRINT_DEBUG
    time_type t1 = get_time();
    if ((tid == 0) && transport.rank() == 0) {
      std::cout << "Time for select A : " << (t1-t0) << std::endl;
    }
#endif

  } // end of operator


  // Active Messages are handled in this function
  void process(const work_item_t& data, int tid) {
    // TODO : simplify message structure after we confirm
    // everything is working. I am suspecting we may need to
    // send a reply after removing vertex from local set
    Vertex u = data.first;
    Vertex v = data.second.first;
    random_t rmap_v = data.second.second.first;
    action_t action = data.second.second.second;

    if (action == REMOVE_VERTEX) {
      if (mis[u] == MIS_UNFIX) {
	if (rmap_v > rmap[u]) {
	  // send a reply -- TODO : We may can simplify the reply
	  work_item_t wi = construct_wi(v, u, rmap[u], REMOVE_REPLY);
	  p_relax_msg->send(wi);
	} else if (rmap_v == rmap[u]) {
	  if (v < u) {
	    // send a reply
	    work_item_t wi = construct_wi(v, u, rmap[u], REMOVE_REPLY);
	    p_relax_msg->send(wi);
	  } else
	    remove_from_local_set(u, tid);
	} else {
	  remove_from_local_set(u, tid);
	}
      }
    } else if (action == REMOVE_REPLY) {
      // just remove the vertex from local set
      remove_from_local_set(u, tid);
    } else if (action == UNFIX_VERTEX) {
      int expected = MIS_UNFIX;
      if (__atomic_compare_exchange_n(&mis[u], &expected, MIS_FIX0, false /*no weak*/, 
				      __ATOMIC_SEQ_CST, __ATOMIC_RELAXED)) {
      }
    } else
      assert(false);
  }

private:
  void add_to_local_set(Vertex v) {
    p_local_mis->push_back(v);
  }

  void remove_from_local_set(Vertex v, int tid) {
    // TODO this could be a map. Not sure which one is going to give
    // better performance. Experiment.
    //boost::mutex::scoped_lock
    //  lock(set_mutex);
    p_deleted_vertices->insert(v, tid);
  }

  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   random_t source_random,
			   action_t action) {
    work_item_t wi(d, std::make_pair(s, std::make_pair(source_random, action)));
    return wi;
  }


  random_t pn; // number of total vertices
  amplusplus::transport& transport;
  Graph& g;
  Bucket_t* p_local_mis;
  OwnerMap* p_owner;
  MISMap& mis;
  RandomMap& rmap;
  DeletedVertices* p_deleted_vertices;
  RelaxMessage* p_relax_msg;
  boost::mutex set_mutex;
  boost::random::mt19937 gen;
  boost::mutex gen_mutex;

};

// SelectAV2 functor generator.
struct select_a_v2_functor_gen {
  template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
  struct select {
    typedef SelectAV2<Graph, MISMap, RandomMap, MessageGenerator> type;
  };
};

//====================== End of SelectA-V2 Functor ==================================//


//======================================Select A Functor - Variation 1 ======================================//
// This variation uses vertex id instead of a random number to find indepenendent sets.
// Functor for Select A - V1 
template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class SelectAV1 {

  // Following definitions must match with the definitions in the graph
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef append_buffer<Vertex> Bucket_t;
  typedef std::pair<Vertex, std::pair<Vertex, std::pair<random_t, action_t> > > work_item_t;
  //  typedef std::set<Vertex> DeletedVertices;
  typedef multi_thread_set<Vertex> DeletedVertices;
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef SelectAV1<Graph, MISMap, RandomMap, MessageGenerator> self_type;

  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function<self_type>, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;


public:
  SelectAV1(Graph& graph, 
	  MISMap& maxmap,
	  RandomMap& random,
	  amplusplus::transport& t) : g(graph), p_local_mis(NULL), 
				      p_owner(NULL),
				      mis(maxmap), rmap(random),
				      p_deleted_vertices(NULL),
				      p_relax_msg(NULL),
				      pn(0), transport(t) {}

  void initialize(Bucket_t* lmis,
		  OwnerMap* owner,
		  DeletedVertices* delv,
		  RelaxMessage* relmsg,
		  unsigned long long numv,
		  int nthreads) {

    p_local_mis = lmis;
    p_owner = owner;
    p_deleted_vertices = delv;
    p_relax_msg = relmsg;
  }

#ifdef MIS_STATS
  void print_stats(int nthreads) {
    /*    auto iterations = rndcalctime[0].size();
    std::cout << "Total number of iterations : " << iterations << std::endl;
    for (int i=0; i < nthreads; ++i) {
      for (int k=0; k < iterations; ++k) {
	std::cout << "Random value calculate Time : " << rndcalctime[i][k] << std::endl;
	std::cout << "Sub-graph compute Time : " << subgcomptime[i][k] << std::endl;
	std::cout << "Independent Set Calculation Time : " << iscalctime[i][k] << std::endl;
	std::cout << "Independent Set Size : " << issize[i][k] << std::endl;
      }
      }*/
  }
#endif


  void operator() (int tid, int nthreads,
		   std::size_t& unfix_set_sz,
		   shared_ptr<amplusplus::detail::barrier>& t_bar) {

    assert(p_local_mis != NULL);
    assert(p_deleted_vertices != NULL);
    assert(p_owner != NULL);
    assert(p_relax_msg != NULL);

#ifdef CRAYPAT
      if (PAT_region_begin ( 3, "av1subset" ) == PAT_API_FAIL) {
	std::cout << "PAT begin failed ! " << std::endl;
	assert(false);
      }
#endif

    BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
      if (get(*p_owner, u) == transport.rank()) {
	if (mis[u] == MIS_UNFIX) {
	  ++unfix_set_sz;
	  rmap[u] = u;
	  add_to_local_set(u);
	}
      }
    }
#ifdef CRAYPAT
    if (PAT_region_end(3) == PAT_API_FAIL) {
      std::cout << "PAT end failed ! " << std::endl;
      assert(false);
    }
#endif


    // synchronize all the nodes
    // thread barrier is impicit ?
    { amplusplus::scoped_epoch epoch(transport); }

#ifdef CRAYPAT
      if (PAT_region_begin ( 4, "av1select" ) == PAT_API_FAIL) {
	std::cout << "PAT begin failed ! " << std::endl;
	assert(false);
      }
#endif


    { 
      amplusplus::scoped_epoch epoch(transport); 

      for (typename Bucket_t::size_type i = tid ; i < p_local_mis->size() ; i += nthreads) {
	Vertex v = (*p_local_mis)[i];
	// There can be many edges between two vertices
	// For algorithm to work we need to make sure
	// each adjacen vertex is counted only once
	// we will have adjacencies set per each thread
	std::set<Vertex> adjacencies;

	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  if (v != u) {
	    if (adjacencies.insert(u).second) {
	      // is u owned by current rank ? then fine, proceed with the algorithm
	      // otherwise send a message to owner with continuation data
	      if (transport.rank() == get(*p_owner, u)) {
		if (mis[u] == MIS_UNFIX) { // This is necessary. Because we are working on the subgraph.
#ifdef PRINT_DEBUG
		  std::cout << "r-v : " << rmap[v] << ", r-u : " << rmap[u] << " local : " 
			    << (transport.rank() == get(*p_owner, u)) << std::endl;
#endif
		  if (rmap[v] > rmap[u]) {
		    // souce must be local to the rank
		    remove_from_local_set(v, tid);
#ifdef PRINT_DEBUG
		    std::cout << "Removing Vertex-a " << v << "- (" << v << "," << u << ")" << std::endl;
#endif
		  } else {
		    remove_from_local_set(u, tid);
#ifdef PRINT_DEBUG
		    std::cout << "Removing Vertex-b " << u << "- (" << v << "," << u << ")" << std::endl;
#endif
		  }
		}
	      } else {
		// u is in a remote rank
		// we need to send v, rmap[v], u, action
		// if rmap[v] > rmap[u], the remote rank needs to send a reply back to current rank : handler
		// if rmap[v] == rmap[u] and v < u remote rank needs to send a reply back to current rank : handler
		// if rmap[v] < rmap[u], the remote rank needs to remove u from its local set, no need to send a reply
		// work_item_t = < target, source, source_random, action >
		work_item_t wi = construct_wi(u, v, rmap[v], REMOVE_VERTEX);
		p_relax_msg->send(wi);
	      }
	    }
	  }
	}
      }
    }

#ifdef CRAYPAT
    if (PAT_region_end(4) == PAT_API_FAIL) {
      std::cout << "PAT end failed ! " << std::endl;
      assert(false);
    }
#endif


  } // end of operator

  // Active Messages are handled in this function
  void process(const work_item_t& data, int tid) {
    // TODO : simplify message structure after we confirm
    // everything is working. I am suspecting we may need to
    // send a reply after removing vertex from local set
    Vertex u = data.first;
    Vertex v = data.second.first;
    random_t rmap_v = data.second.second.first;
    action_t action = data.second.second.second;

    if (action == REMOVE_VERTEX) {
      if (mis[u] == MIS_UNFIX) {
	if (rmap_v > rmap[u]) {
	  // send a reply -- TODO : We may can simplify the reply
	  work_item_t wi = construct_wi(v, u, rmap[u], REMOVE_REPLY);
	  p_relax_msg->send(wi);
	} else {
	  remove_from_local_set(u, tid);
	}
      }
    } else if (action == REMOVE_REPLY) {
      // just remove the vertex from local set
      remove_from_local_set(u, tid);
    } else if (action == UNFIX_VERTEX) {
      int expected = MIS_UNFIX;
      if (__atomic_compare_exchange_n(&mis[u], &expected, MIS_FIX0, false /*no weak*/, 
				      __ATOMIC_SEQ_CST, __ATOMIC_RELAXED)) {
      }
    } else
      assert(false);
  }

private:
  void add_to_local_set(Vertex v) {
    p_local_mis->push_back(v);
  }

  void remove_from_local_set(Vertex v, int tid) {
    // TODO this could be a map. Not sure which one is going to give
    // better performance. Experiment.
    //boost::mutex::scoped_lock
    //  lock(set_mutex);
    p_deleted_vertices->insert(v, tid);
  }

  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   random_t source_random,
			   action_t action) {
    work_item_t wi(d, std::make_pair(s, std::make_pair(source_random, action)));
    return wi;
  }


  random_t pn; // number of total vertices
  amplusplus::transport& transport;
  Graph& g;
  Bucket_t* p_local_mis;
  OwnerMap* p_owner;
  MISMap& mis;
  RandomMap& rmap;
  DeletedVertices* p_deleted_vertices;
  RelaxMessage* p_relax_msg;
  boost::mutex set_mutex;
};


// SelectAV1 functor generator.
struct select_a_vertex_functor_gen {
  template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
  struct select {
    typedef SelectAV1<Graph, MISMap, RandomMap, MessageGenerator> type;
  };
};

//====================== End of SelectA-V1 Functor ==================================//


//======================================Select B Functor ======================================//
//boost::mutex inc_mutex;
// This is the Luby's B algorithm
// Functor for Select B
template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class SelectB {

  // Following definitions must match with the definitions in the graph
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef append_buffer<Vertex> Bucket_t;
  typedef std::pair<Vertex, std::pair<Vertex, std::pair<random_t, action_t> > > work_item_t;
  //  typedef std::set<Vertex> DeletedVertices;
  typedef multi_thread_set<Vertex> DeletedVertices;
  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef SelectB<Graph, MISMap, RandomMap, MessageGenerator> self_type;

  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function<self_type>, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;

  typedef unsigned long degree_t;
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef iterator_property_map<typename std::vector<degree_t>::iterator, VertexIndexMap>  DegreeMap;



public:
  SelectB(Graph& graph, 
	  MISMap& maxmap,
	  RandomMap& random,
	  amplusplus::transport& t) : g(graph), p_local_mis(NULL), 
				      p_owner(NULL),
				      mis(maxmap), rmap(random),
				      p_deleted_vertices(NULL),
				      p_relax_msg(NULL),
				      pn(0), transport(t) {}

  void initialize(Bucket_t* lmis,
		  OwnerMap* owner,
		  DeletedVertices* delv,
		  RelaxMessage* relmsg,
		  unsigned long long numv,
		  int nthreads) {

    p_local_mis = lmis;
    p_owner = owner;
    p_deleted_vertices = delv;
    p_relax_msg = relmsg;

    p_degrees = new std::vector<degree_t>(num_vertices(g), 0);
    p_degree_map = new DegreeMap(p_degrees->begin(), get(vertex_index, g));

#ifdef MIS_STATS
    //    selecttime.resize(nthreads);
    //mergetime.resize(nthreads);
    rndcalctime.resize(nthreads);
    subgcomptime.resize(nthreads);
    iscalctime.resize(nthreads);
    localsetsz.resize(nthreads);
    removesetsz.resize(nthreads);
#endif
    

  }

#ifdef MIS_STATS
  void print_stats(int nthreads) {
    auto iterations = rndcalctime[0].size();
    std::cout << "Total number of iterations : " << iterations << std::endl;
    for (int i=0; i < nthreads; ++i) {
      for (int k=0; k < iterations; ++k) {
	//std::cout << "tid : " << i << ", ite : " << k << ", Random value calculate Time : " << rndcalctime[i][k] << std::endl;
	//std::cout << "tid : " << i << ", ite : " << k << ", Sub-graph compute Time : " << subgcomptime[i][k] << std::endl;
	std::cout << "tid : " << i << ", ite : " << k << ", Independent Set Calculation Time : " << iscalctime[i][k] << std::endl;
	std::cout << "tid : " << i << ", ite : " << k << ", Local Independent Set Size : " << localsetsz[i][k] << std::endl;
	std::cout << "tid : " << i << ", ite : " << k << ", Remove Set Size : " << removesetsz[i][k] << std::endl;
      }
    }

    std::cout << "========================== Aggregated Independent Set Sizes ==============================" << std::endl;
    //aggregate over threads
    std::vector<unsigned long> aggrlocalset;
    std::vector<unsigned long> aggrremoveset;
    aggrlocalset.resize(iterations, 0);
    aggrremoveset.resize(iterations, 0);

    for (int k=0; k < iterations; ++k) {
      for (int i=0; i < nthreads; ++i) {
	aggrlocalset[k] += localsetsz[i][k];
	aggrremoveset[k] += removesetsz[i][k];
      } 
    }

    for (int i=0; i < iterations; ++i) {
      std::cout << "ite : " << i << ", Assumed Independent Set Size : " << aggrlocalset[i] << std::endl;
      std::cout << "ite : " << i << ", Deleted Set Size : " << aggrremoveset[i] << std::endl;
    }


  }
#endif



  ~SelectB() {
    // de-allocate iterator property map
    // and the vector
    delete p_degrees;
    p_degrees = NULL;

    delete p_degree_map;
    p_degree_map = NULL;
  }

  int coin(Vertex i) {
    degree_t degree = (*p_degree_map)[i];
    assert(degree >= 0);

    if (degree == 0)
      return 1;
 
    double probability = (double)(1/(double)(2*degree));

    #ifdef PRINT_DEBUG
    std::cout << "Vertex : " << i << ", degree : " 
	      << degree << ", probability : " 
	      << probability << std::endl;
    #endif

    double probabilities[] = {
      (1-probability), probability
    };

    boost::random::discrete_distribution<> dist(probabilities);
    int ret = 0;
    {
      boost::mutex::scoped_lock
	lock(rng_mutex);
      ret = dist(gen); 
    }

    return ret;
  }


  void operator() (int tid, int nthreads,
		   std::size_t& unfix_set_sz,
		   shared_ptr<amplusplus::detail::barrier>& t_bar) {

    assert(p_local_mis != NULL);
    assert(p_deleted_vertices != NULL);
    assert(p_owner != NULL);
    assert(p_relax_msg != NULL);
    assert(p_degree_map != NULL);

    if (tid == 0)
      setx.clear();

    t_bar->wait();

    // There can be many edges between two vertices
    // For algorithm to work we need to make sure
    // each adjacen vertex is counted only once
    // we will have adjacencies set per each thread
    std::set<Vertex> adjacencies;
    // calculate degrees relative to subgraph
#ifdef MIS_STATS
    time_type asubst = get_time();
#endif

    { 
      amplusplus::scoped_epoch epoch(transport); 

      BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
	if (get(*p_owner, u) == transport.rank()) {
	  if (mis[u] == MIS_UNFIX) {
	    //increment the unfix vertices
	    // This is to check whether subgraph is empty
	    // in the main algorithm
	    ++unfix_set_sz;

	    (*p_degree_map)[u] = 0;
	    //	    BGL_FORALL_ADJ_T(u, v, g, Graph) {
	    BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
	      Vertex v = target(e, g);
	      // ignore self-loops
	      if (v != u) {
		// count adjacent vertex only if it is not counted
		// already
		if (adjacencies.insert(v).second) {// insert adjacen vertex to adjacencies set

		  if (get(*p_owner, v) == transport.rank()) {
		    if (mis[v] == MIS_UNFIX) {
		      //		      boost::unique_lock<boost::mutex> scoped_lock(inc_mutex);
		      ++(*p_degree_map)[u];
		    }
		  } else {
		    // send a message to remote vertex to check
		    // whether vertex is in MIS or not
		    // if it is not in mis
		    work_item_t wi = construct_wi(v /*destination*/, 
						  u /*source*/, 
						  0, 
						  IN_MIS_QUERY);
		    p_relax_msg->send(wi);

		  }
		}
	      }
	    }
	  }
	}
	// clear adjacent vertices for current source vertex
	adjacencies.clear();
      }
    }

#ifdef MIS_STATS
    time_type esubst = get_time();
    subgcomptime[tid].push_back(esubst - asubst);
#endif

#ifdef MIS_STATS
  time_type srndtm = get_time();
#endif

    // by the end of this epoch we have calculated  all  the
    // vertex degree for vertices that are not decided whether
    // in mis or not in mis (MIS_UNFIX).
    BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
      if (get(*p_owner, u) == transport.rank()) {
	if (mis[u] == MIS_UNFIX) {
	  if (coin(u) == 1) {
	    add_to_local_set(u);
	  }
	}
      }
    }

#ifdef MIS_STATS
  time_type erndtm = get_time();
  rndcalctime[tid].push_back(erndtm-srndtm);
  localsetsz[tid].push_back(p_local_mis->size());
#endif

    // synchronize all the nodes
    // thread barrier is impicit ?
    { amplusplus::scoped_epoch epoch(transport); }
#ifdef MIS_STATS
    time_type t0 = get_time();
#endif

    { 
      amplusplus::scoped_epoch epoch(transport); 

      for (typename Bucket_t::size_type i = tid ; i < p_local_mis->size() ; i += nthreads) {
	Vertex v = (*p_local_mis)[i];
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  if (v != u) {
	
	    // is u owned by current rank ? then fine, proceed with the algorithm
	    // otherwise send a message to owner with continuation data
	    if (transport.rank() == get(*p_owner, u)) {
	      if (is_in_setx(u)) { // This is necessary. Because we are working on the subgraph.
		// we know that v is in set X
		// but we do not know whether u is in set X
		// searching in an append buffer is inefficient. Therefore, 
		// we calculate coin value again for u. if coin(u) == 1 we know
		// u is also in X. In that case, based on the degree we 
		// need to remove a vertex

		if ((*p_degree_map)[v] <= (*p_degree_map)[u]) {
                  remove_from_local_set(v, tid);
		} else { 
                  remove_from_local_set(u, tid);
		}
	      }
	    } else {
	      // u is in a remote rank
	      // we need to send v, (*p_degree_map)[v], u, action
	      // if (*p_degree_map)[v] <= (*p_degree_map)[u], the remote rank needs to send a reply back to current rank : handler
	      // if (*p_degree_map)[v] > (*p_degree_map)[u], the remote rank needs to remove u from its local set, no need to send a reply
	      work_item_t wi = construct_wi(u, v, (*p_degree_map)[v], REMOVE_VERTEX);
	      p_relax_msg->send(wi);
	    }
	  }
	}
      }
    }

#ifdef MIS_STATS
    time_type t1 = get_time();
    iscalctime[tid].push_back(t1-t0);
    removesetsz[tid].push_back(p_deleted_vertices->size());
#endif


  } // end of operator

  // Active Messages are handled in this function
  void process(const work_item_t& data, int tid) {

    // TODO : simplify message structure after we confirm
    // everything is working. I am suspecting we may need to
    // send a reply after removing vertex from local set
    Vertex u = data.first;
    Vertex v = data.second.first;
    state_t degree_v = data.second.second.first;
    action_t action = data.second.second.second;

    if (action == REMOVE_VERTEX) {
      if (is_in_setx(u)) {
	// is u in X ?
	if (degree_v <= (*p_degree_map)[u]) {
	  // In this case v must be removed, therefore send a reply back to v
	  work_item_t wi = construct_wi(v, u, (*p_degree_map)[u], REMOVE_REPLY);
	  p_relax_msg->send(wi);
	} else {
          remove_from_local_set(u, tid);
	}
      }
    } else if (action == REMOVE_REPLY) {
      // just remove the vertex from local set
      remove_from_local_set(u, tid);
    } else if (action == UNFIX_VERTEX) {
      int expected = MIS_UNFIX;
      if (__atomic_compare_exchange_n(&mis[u], &expected, MIS_FIX0, false /*no weak*/, 
				      __ATOMIC_SEQ_CST, __ATOMIC_RELAXED)) {
      }

    } else if (action == IN_MIS_QUERY) { 
      // NOTE : If it is decided that u is in MIS
      // or not we will not send a reply. 
      if (mis[u] == MIS_UNFIX) {
	work_item_t wi = construct_wi(v, u, 1, IN_MIS_REPLY);
	p_relax_msg->send(wi);
      }
    } else if (action == IN_MIS_REPLY) {
      // increase the degree
      // following statement should be executed in a lock
      //      boost::unique_lock<boost::mutex> scoped_lock(inc_mutex);
      ++(*p_degree_map)[u];
    } else
      assert(false);
  }

private:
  void add_to_local_set(Vertex v) {
    p_local_mis->push_back(v);
    
    // we also need to push it to set cos
    // searching in append buffer is inefficient
    {
      boost::mutex::scoped_lock
	lock(set_mutex);
      assert(mis[v] == MIS_UNFIX);
      setx.insert(v);
    }    
  }

  // returns true if v is in setx
  bool is_in_setx(Vertex v) {
#ifdef PRINT_DEBUG
    if (setx.find(v) != setx.end())
      assert(mis[v] == MIS_UNFIX);
#endif
    return (setx.find(v) != setx.end());
  }

  void remove_from_local_set(Vertex v, int tid) {
    // TODO this could be a map. Not sure which one is going to give
    // better performance. Experiment.
    //boost::mutex::scoped_lock
    //  lock(set_mutex);
    p_deleted_vertices->insert(v, tid);
  }

  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   random_t source_random,
			   action_t action) {
    work_item_t wi(d, std::make_pair(s, std::make_pair(source_random, action)));
    return wi;
  }


  random_t pn; // number of total vertices
  amplusplus::transport& transport;
  Graph& g;
  Bucket_t* p_local_mis;
  OwnerMap* p_owner;
  MISMap& mis;
  std::vector<degree_t>* p_degrees;
  std::set<Vertex> setx;
  DegreeMap* p_degree_map;
  RandomMap& rmap;
  DeletedVertices* p_deleted_vertices;
  RelaxMessage* p_relax_msg;
  boost::mutex set_mutex;
  boost::mutex rng_mutex;
  boost::mt19937 gen;
#ifdef MIS_STATS
  std::vector<std::vector<time_t> > rndcalctime;
  std::vector<std::vector<time_t> > mergetime;
  std::vector<std::vector<time_t> > subgcomptime;
  std::vector<std::vector<time_t> > iscalctime;
  std::vector<std::vector<unsigned long> > localsetsz;
  std::vector<std::vector<unsigned long> > removesetsz;
#endif

};


// SelectB functor generator.
struct select_b_functor_gen {
  template<typename Graph, 
	 typename MISMap,
	 typename RandomMap,
	 typename MessageGenerator =
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
  struct select {
    typedef SelectB<Graph, MISMap, RandomMap, MessageGenerator> type;
  };
};

//====================== End of SelectB Functor ==================================//


template<typename Graph, typename MISMap, typename RandomMap, 
	 typename SelectFunctorGen=boost::graph::distributed::select_a_functor_gen,
	 typename MessageGenerator = 
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class luby_mis {

  //  typedef luby_mis<Graph, MISMap, MessageGenerator> 
  //  self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;
  typedef typename std::map<Vertex, std::set<Vertex> > VertexSetType;
  typedef append_buffer<Vertex> Bucket_t;

  // work_item_t = < target, source, source_random, action >
  typedef std::pair<Vertex, std::pair<Vertex, std::pair<random_t, action_t> > > work_item_t;

  //  typedef std::set<Vertex> DeletedVertices;
  typedef multi_thread_set<Vertex> DeletedVertices;

  // The SelectFunctor type
  typedef typename SelectFunctorGen::template select<Graph, MISMap, RandomMap, MessageGenerator>::type SelectFunctor;

  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function<SelectFunctor>, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;
  
public:
  luby_mis(Graph& g,
	   MISMap& mismap,
	   RandomMap& random,
	   unsigned long long numv, // number of vertices
	   amplusplus::transport &t,
	   int offset,
	   MessageGenerator message_gen = 
	   MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
      g(g), transport(t), mis(mismap), rmap(random),
      n(numv), owner(get(vertex_owner, g)), 
      relax_msg(message_gen, transport, owner_from_pair<OwnerMap, work_item_t>(owner),
		amplusplus::idempotent_combination(minimum_pair_first())),
    core_offset(offset),
    select_function(g, mismap, random, t),
    deleted_vertices(t.get_nthreads())
  {
    initialize();
  }

  void operator() (int tid) { run(tid); }
  void run(int tid = 0);
  time_type get_start_time() {
    return start_time;
  }

#ifdef MIS_STATS
  void print_stats() {
    select_function.print_stats(transport.get_nthreads());
  }
#endif



private:
  void initialize();
  void process(const work_item_t& data, int tid);

  // TODO remove
  work_item_t construct_wi(Vertex d,
			   Vertex s,
			   random_t source_random,
			   action_t action) {
    work_item_t wi(d, std::make_pair(s, std::make_pair(source_random, action)));
    return wi;
  }

  // TODO remove
  /*void add_to_local_set(Vertex v) {
    local_mis.push_back(v);
  }

  // TODO remove
  void remove_from_local_set(Vertex v) {
    // TODO this could be a map. Not sure which one is going to give
    // better performance. Experiment.
    boost::mutex::scoped_lock
      lock(set_mutex);
    deleted_vertices.insert(v);
    }*/

  void clear_containers() {
    local_mis.clear();
    deleted_vertices.clear();

    // clear ghost cells
    mis.clear();
    rmap.clear();
  }

  //void select(int tid, int nthreads);

  const int dummy_first_member_for_init_order;
  const Graph& g;
  amplusplus::transport& transport;
  MISMap& mis;
  RandomMap& rmap;
  unsigned long long n; // number of total vertices
  OwnerMap owner;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  RelaxMessage relax_msg;
  time_type start_time;
  SelectFunctor select_function;
  DeletedVertices deleted_vertices;
  boost::mutex set_mutex;
  Bucket_t local_mis;
  int core_offset;

};

#define MIS_PARMS                                   \
      typename Graph, typename MISMap, typename RandomMap, typename SelectFunctorGen, typename MessageGenerator

#define MIS_TYPE                                    \
      luby_mis<Graph, MISMap, RandomMap, SelectFunctorGen, MessageGenerator>


template<MIS_PARMS>
void
MIS_TYPE::initialize() {
  // Set the handler
  relax_msg.set_handler(processing_function<SelectFunctor>(&select_function));

  // initialize select_function
  select_function.initialize(&local_mis,
			     &owner,
			     &deleted_vertices,
			     &relax_msg,
			     n,
			     transport.get_nthreads());

  BGL_FORALL_VERTICES_T(u, g, Graph) {
    if (get(owner, u) == transport.rank()) {
      mis[u] = MIS_UNFIX;
    }
  }
}



template<MIS_PARMS>
void
MIS_TYPE::run(int tid) {
  int ite_count = 0;
  AMPLUSPLUS_WITH_THREAD_ID(tid) {

#ifdef CRAYPAT
      if (PAT_region_begin ( 1, "lubyrun" ) == PAT_API_FAIL) {
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

    std::size_t all_sizes = 0;
    while(true) {
      std::size_t unfix_set_sz = 0;
#ifdef PRINT_DEBUG
      time_type t1 = get_time();
#endif

      //#ifdef CRAYPAT
      //if (PAT_region_begin ( 2, "lubyselect" ) == PAT_API_FAIL) {
      //	std::cout << "PAT begin failed ! " << std::endl;
      //	assert(false);
      //}
      //#endif


      select_function(tid, nthreads, unfix_set_sz, t_bar);

      //#ifdef CRAYPAT
	//if (PAT_region_end(2) == PAT_API_FAIL) {
      //std::cout << "PAT end failed ! " << std::endl;
      //assert(false);
      //}
      //#endif


      t_bar->wait();

#ifdef PRINT_DEBUG
      time_type t2 = get_time();
#endif

      {
	amplusplus::scoped_epoch_value epoch(transport, unfix_set_sz, all_sizes);

#ifdef PRINT_DEBUG
	if (tid == 0 && (transport.rank() == 0)) {
	  std::cout << "local mis size : " << local_mis.size() << std::endl;
	  std::cout << "deleted vertices size : " << deleted_vertices.size() << std::endl;
	  std::cout << "Independent set size : " << (local_mis.size() - deleted_vertices.size()) << std::endl;
	}

	if ((tid == 0) &&
	    (transport.rank() == 0)) {
	  std::cout << "Unfix set size : " << all_sizes
		    << ", select time : " << (t2-t1) << std::endl;
	}
#endif

	for (typename Bucket_t::size_type i = tid ; i < local_mis.size() ; i += nthreads) {
	  Vertex v = local_mis[i];
	  // check whether v is deleted
	  if (!deleted_vertices.exists(v, tid)) {
	    // vertex v is not deleted
	    assert(mis[v] == MIS_UNFIX);
	    mis[v] = MIS_FIX1;
	    BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	      Vertex u = target(e, g);
	      if (v != u) { // ignore self loops
		work_item_t wi = construct_wi(u, v, rmap[v], UNFIX_VERTEX);
		relax_msg.send(wi);
	      }
	    }
	  }
	}

#ifdef USE_PARALLEL_FORALL_EDGES_MACRO
	if (tid == 0) { // TODO parallelize this ! -- probably we need to reply set with a map
	  typename std::set<Vertex>::iterator ite = local_mis.begin();
	  for (; ite != local_mis.end(); ++ite) {
	    assert(mis[(*ite)] == MIS_UNFIX);

	    // add vertices in the set to mis
	    mis[(*ite)] = MIS_FIX1;
	    BGL_FORALL_OUTEDGES_T((*ite), e, g, Graph) {
	      Vertex u = target(e, g);
	      work_item_t wi = construct_wi(u, (*ite), UNFIX_VERTEX, tid);
	      relax_msg.send(wi);
	    }
	  }
	}
#endif
      }

      if (all_sizes == 0)
	break;

      if (tid == 0) {
	clear_containers(); // clear all local independent sets before proceeding further
	++ite_count;
      }
      
      // wait till everyone clear local mis sets
      { amplusplus::scoped_epoch epoch(transport); }
    }

#ifdef CRAYPAT
    if (PAT_region_end(1) == PAT_API_FAIL) {
      std::cout << "PAT end failed ! " << std::endl;
      assert(false);
    }
#endif


  }

  if ((tid == 0) && transport.rank() == 0)
    std::cout << "Number of iterations : " << ite_count << std::endl;

}



template<typename SelectFunctor>
struct
processing_function {
  
  processing_function() : select_fn(NULL) {}
  processing_function(SelectFunctor* sf) : select_fn(sf) {}
  
  template<typename work_item_t>
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
#ifdef PRINT_DEBUG
    std::cout << "Handler called in tid : " << tid << std::endl;
#endif
    select_fn->process(data, tid);
  }

protected:
  SelectFunctor* select_fn;
};

}}}
#endif
