// Copyright (C) 2017 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Triangle Counting Algortihms================//
// Driver for BFS Family Algorithms
//===========================================================//
#ifndef BFS_EXECUTOR
#define BFS_EXECUTOR

#include <iostream>
#include "../common/synthetic_generator.hpp"
#include "../common/parser.hpp"
#include "../common/executor.hpp"
#include "../common/vertex_permutations.hpp"
#include "../common/spinlock.hpp"
#include <boost/graph/distributed/bfs_chaotic.hpp>
#include <boost/random/linear_congruential.hpp>
#include "graph_stats.hpp"


enum bfs_algorithm {
  chaotic,
  level_sync,
  k_level,
  dc,
  stats
};

class bfs_instance_params {

public:
  bfs_algorithm algorithm;
  bool verify;
  uint64_t source = std::numeric_limits<uint64_t>::max();

  bfs_instance_params(bfs_algorithm alg, bool v, int level, uint64_t s) : algorithm(alg), 
									  verify(v),
									  source(s){}

  bfs_instance_params(bfs_algorithm alg, bool v, int level) : algorithm(alg), 
							      verify(v){}
							
  bfs_instance_params(bfs_algorithm alg, bool v) : algorithm(alg), 
						   verify(v){}
					
							


  void print() {
    std::cout << "Algorithm : " << get_algorithm() << std::endl;
    std::cout  << "Source : " << source << std::endl;
  }

  std::string get_algorithm() {
    if (algorithm == chaotic) 
      return "BFS Chaotic";
    else if(algorithm == level_sync)
      return "BFS Level Sync";
    else if(algorithm == k_level)
      return "BFS k-level";
    else if(algorithm == dc)
      return "BFS DC";
    else if(algorithm == stats)
      return "stats";
    else
      return "Invalid Algorithm!";
  }

};

class bfs_params {

private:
  std::vector<bfs_instance_params> instance_params;
  std::vector<bfs_algorithm> algorithms;
  std::vector<uint64_t> sources;
  bool verify = false;

public:
  bool parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      bfs_algorithm algorithm;
      std::string arg = argv[i];
      if (arg == "--run-bfs-chaotic") {
	algorithm = chaotic;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run-bfs-level") {
	algorithm = level_sync;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run-bfs-klevel") {
	algorithm = k_level;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run-bfs-dc") {
	algorithm = dc;
	algorithms.push_back(algorithm);
      }

      if (arg == "--stats") {
	algorithm = stats;
	algorithms.push_back(algorithm);
      }

      if (arg == "--sources") {
	sources = extract_params<uint64_t>( argv[i+1] );
      }

      if (arg == "--verify") {
	verify = true;
      }
    }
  }

  void print() {
    std::cout << "[INFO] Algorithms : [";
    BOOST_FOREACH(bfs_algorithm algo, algorithms) {
      std::cout << get_algorithm(algo) << ", ";
    }
    std::cout << "]";

    std::cout << "[INFO] Verify : " << verify << std::endl;
  }

  std::string get_algorithm(bfs_algorithm algorithm) {
    if (algorithm == chaotic) 
      return "BFS Chaotic";
    else if(algorithm == level_sync)
      return "BFS Level Sync";
    else if(algorithm == k_level)
      return "BFS k-level";
    else if(algorithm == dc)
      return "BFS DC";
    else if(algorithm == stats)
      return "Graph Stats";
    else
      return "Invalid Algorithm!";
  }


  std::vector<bfs_instance_params>&
  get_instance_params() {   
    if (sources.size() == 0) {
      BOOST_FOREACH(bfs_algorithm alg, algorithms) {
	instance_params.push_back(bfs_instance_params(alg, verify));
      }
    } else {
      BOOST_FOREACH(uint64_t s, sources) {
	BOOST_FOREACH(bfs_algorithm alg, algorithms) {
	  instance_params.push_back(bfs_instance_params(alg, verify, std::numeric_limits<int>::max(), s));
	}
      }
    }

    return instance_params;
  }

};


class BFSExecutor {

private:
  int abort_level;
  uint64_t useful_work;
  uint64_t rejected_work;
  uint64_t invalidated_work;

  template <typename Graph, typename Vertex, typename graph_create_params>
  Vertex  get_random_source(Graph& g, graph_create_params& gparams) {
    boost::uniform_int<Vertex> rand_vertex(0, gparams.n-1);
    boost::random::rand48 synch_gen;
    synch_gen.seed(static_cast<unsigned int>(std::time(0)));

    return static_cast<Vertex>(boost::vertex(rand_vertex(synch_gen), g));

  }

  template <typename Graph, typename DistanceMap>
  bool verify(Graph& g, 
	      DistanceMap& distance) {

    distance.set_consistency_model(boost::parallel::cm_forward);
    distance.set_max_ghost_cells(0);

    if (g.transport().rank() == 0) std::cout<<"Verifying results......";

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
	if (get(distance, target(e, g)) > boost::closed_plus<Distance>()(get(distance, source(e, g)), 1))
	  std::cout << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g)) << "->"
		    << get(get(vertex_local, g), target(e, g)) << "@" << get(get(vertex_owner, g), target(e, g)) << "  weight = "
		    << get(weight, e)
		    << "  distance(" << get(get(vertex_local, g), source(e, g)) << "@" << get(get(vertex_owner, g), source(e, g))
		    << ") = " << get(distance, source(e, g)) << "  distance(" << get(get(vertex_local, g), target(e, g)) << "@" 
		    << get(get(vertex_owner, g), target(e, g)) << ") = " << get(distance, target(e, g)) << std::endl;
#else
	if(get(distance, target(e, g)) > boost::closed_plus<Distance>()(get(distance, v), 1)) std::abort();
#endif
      }
    }
    if (g.transport().rank() == 0) std::cout << "Verified." << std::endl;
    distance.clear(); // Clear memory used by ghost cells

    return true;
  }


  template <typename Graph, 
	    typename MessageGenerator,
	    typename graph_create_params>
  time_type
  run_bfs_chaotic(amplusplus::transport& trans, 
		  amplusplus::transport& barrier_trans, 
		  Graph& g,  
		  MessageGenerator msg_gen, 
		  graph_create_params& gparams,
		  instance_params& runtime_params,
		  bfs_instance_params& bfs_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    typedef boost::iterator_property_map<typename std::vector<Distance>::iterator, VertexIndexMap>  DistanceMap;

    std::vector<Distance> distvec(num_vertices(g), std::numeric_limits<Distance>::max());
    DistanceMap dist_map(distvec.begin(), get(boost::vertex_index, g));

    if (trans.rank() == 0)
      std::cout << "Initializing BFS chaotic ..." << std::endl;

    boost::graph::distributed::bfs_chaotic<Graph, DistanceMap, MessageGenerator>
      D(g, dist_map, trans, sched_getcpu(), msg_gen);

    // globale edges traversed
    uint64_t global = 0;

    time_type exectime = 0;

    Vertex current_source;

    while(true) {

      global = 0;
      exectime = 0;

      trans.set_nthreads(1);

      { amplusplus::scoped_epoch epoch(barrier_trans); }

      if (bfs_params.source != std::numeric_limits<uint64_t>::max())
	current_source = bfs_params.source;
      else
	current_source = get_random_source<Graph, Vertex, graph_create_params>(g, gparams);

      // Many threads now
      trans.set_nthreads(runtime_params.threads);

      if (trans.rank() == 0)
	std::cout << "Invoking algorithm ..." << std::endl;

      boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
      for (int i = 0; i < runtime_params.threads - 1; ++i) {
	boost::thread thr(boost::ref(D), i + 1);
	threads[i].swap(thr);
      }
	  
      D.set_source(current_source);
      D.set_abort_level(abort_level);
      D(0);
    
      for (int i = 0; i < runtime_params.threads - 1; ++i)
	threads[i].join();
	  
      time_type end = get_time();

      if (trans.rank() == 0)
	std::cout << "Algorithm done ..." << std::endl;

      time_type start = D.get_start_time();

      exectime = (end-start);

      // Back to one thread
      trans.set_nthreads(1);
      clear_thread_core_data();

      std::cout << "verify" << std::endl;
      // verification
      if (bfs_params.verify) {
	std::cout << "inside" << std::endl;
	assert(verify(g, dist_map));
      }

      // local edges traversed
      uint64_t local = D.get_traversed_edges();
    
      MPI_Reduce(&local, &global, 
		 1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);

      std::cout << "Number of edges traversed : " << global << std::endl;
      // if algorithm did not traverse at least half of the edges 
      // continue with an another source
      /*if ((global / 8) > gparams.n) {
	std::cout << "[INFO] Algorithm did not traverse at least half of the vertices. Therefore, continuing with another source " << std::endl;

	std::cout << "source : " << bfs_params.source << std::endl;
	if (abort_level == std::numeric_limits<uint64_t>::max())
	  continue;
	  }*/

      break;
    }
    if (_RANK==0) {
      std::cout << "[INFO] Source " << current_source << " Total edges traversed : " << global << ", time : " << exectime
		<< ", TEPS : " << (double)(global/exectime) << std::endl;
    }

#ifdef COLLECT_STATS
    std::cout << "Level : " << abort_level << std::endl;
    D.print_stats();
    work_stats_t localstats = D.get_local_work_stats();
    // useful, invalidated, rejected
    uint64_t luseful = std::get<0>(localstats);
    uint64_t linvalidated = std::get<1>(localstats);
    uint64_t lrejected = std::get<2>(localstats);

    std::cout << "luseful=" << luseful << std::endl;

    MPI_Reduce(&luseful, &useful_work, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&linvalidated, &invalidated_work, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&lrejected, &rejected_work, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
#endif


    return exectime;

  }

public:
  BFSExecutor():abort_level(std::numeric_limits<int>::max()), useful_work(0),
		rejected_work(0), invalidated_work(0){}

  void set_abort_level(int alevel) {
    abort_level = alevel;
  }

  uint64_t get_useful_work() {
    return useful_work;
  }

  uint64_t get_rejected_work() {
    return rejected_work;
  }

  uint64_t get_invalidated_work() {
    return invalidated_work;
  }

  void reset_work_stats() {
    abort_level = std::numeric_limits<int>::max();
    useful_work = 0;
    rejected_work = 0;
    invalidated_work = 0;
  }

  template <typename Graph, typename MessageGenerator, typename graph_create_params>
  time_type operator()(const Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       bfs_instance_params& bfs_params) {

    amplusplus::transport barrier_trans = trans.clone();
    if (bfs_params.algorithm == chaotic) {
      return run_bfs_chaotic(trans, barrier_trans,
			     g, msg_gen, gparams,
			     runtime_params,
			     bfs_params);
    } else if (bfs_params.algorithm == stats) {
      graph_stats<Graph> gstats;
      time_type t1 = get_time();
      gstats.print_degree_distribution(g);
      time_type t2 = get_time();
      return (t2-t1);
    } else {
      std::cout << "[ERROR] Only BFS chaotic is supported." << std::endl;
      assert(false);
    }

  }
};
#endif
