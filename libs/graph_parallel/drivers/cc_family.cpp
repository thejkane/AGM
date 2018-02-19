// Copyright (C) 2017 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== CC Algortihms================//
// Driver for CCfamily of algorithms.
//===========================================================//


#include <iostream>
#include "common/synthetic_generator.hpp"
#include "common/parser.hpp"
#include "common/executor.hpp"
#include "common/vertex_permutations.hpp"
#include "common/ampp_perf_counters.hpp"
#include "common/work_stats.hpp"

#include <boost/graph/distributed/cc_dc.hpp>
#include <boost/graph/distributed/cc_chaotic.hpp>
#include <boost/graph/distributed/level_sync_cc.hpp>
#include <boost/graph/distributed/delta_stepping_cc.hpp>

#define NODE_PRIORITY_Q_GEN boost::graph::distributed::node_priority_queue_gen
#define NUMA_PRIORITY_Q_GEN boost::graph::distributed::numa_priority_queue_gen
#define THREAD_PRIORITY_Q_GEN boost::graph::distributed::cc_default_priority_queue_gen

enum cc_agm {
  delta,
  level,
  chaotic,
  undefined_agm
};

enum cc_eagm {
  global,
  node,
  numa,
  thread,
  undefined_eagm
};

// Shiloach-Vishkin variations
// TODO
enum cc_sv {
  undefined_sv
};



class cc_instance_params {
public:
  cc_agm agm;
  cc_eagm eagm;
  cc_sv sv;
  id_distribution_t id_distribution;
  bool verify;
  int delta_val;
  agm_work_stats work_stats;
  
  cc_instance_params(int threads,
                     cc_agm _agm,
		     cc_eagm _eagm,
		     cc_sv _sv,
		     id_distribution_t& idd,
		     bool v,
		     int d):
    agm(_agm), 
    eagm(_eagm), 
    sv(_sv),
    id_distribution(idd), 
    verify(v),
    delta_val(d),
    work_stats(threads){}

  ~cc_instance_params() {
  }

  std::string get_agm() {
    if (agm == delta)
      return "Delta";
    else if (agm == level)
      return "Level";
    else if (agm == chaotic)
      return "Chaotic";
    else
      return "Invalid Algorithm";
  }

  std::string get_eagm() {
    if (eagm == global)
      return "Global";
    else if (eagm == node)
      return "Node";
    else if (eagm == numa)
      return "NUMA";
    else if (eagm == thread)
      return "Thread";
    else
      return "Invalid EAGM";
  }

  std::string get_sv() {
    //TODO
    return "";
  }

  void print() {
    std::cout << "Algorithm : " << get_algorithm() << std::endl;

    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "delta val : " << delta_val << std::endl;
    std::cout << "verify : " << verify << std::endl;
  }

  std::string get_algorithm() {
    std::string algorithm = "[";

    if (agm != undefined_agm) {
      algorithm += get_agm();
      algorithm += ", ";
      algorithm += get_eagm();
    } else {
      algorithm += get_sv();
    }
    algorithm += "]";

    return algorithm;
  }

  void print_summary() {
    work_stats.print_summary();
  }

};


class cc_params {

private:
  std::vector<cc_instance_params> instance_params;
  std::vector<cc_agm> agms;
  std::vector<cc_eagm> eagms;
  // TODO std::vector<cc_sv> svs;
  id_distribution_t id_distribution = horizontal; // default
  bool verify = false;
  std::vector<int> deltas;
  int threads;

public:
  bool parse(int argc, char* argv[]){
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--agm") {
	if (strcmp(argv[i+1],"delta") == 0) {
	  agms.push_back(delta);
	} else if (strcmp(argv[i+1],"level") == 0) {
	  agms.push_back(level);
	} else if (strcmp(argv[i+1],"chaotic") == 0) {
	  agms.push_back(chaotic);
	} else {
	  std::cout << "[ERROR] Invalid AGM. Available : delta, choatic, level" << std::endl;
	  return false;
	}
      }

      if (arg == "--delta-val") {
	deltas = extract_params<int>( argv[i+1] );
      }

      if (arg == "--threads") {
        threads = boost::lexical_cast<int>( argv[i+1] );
      }

      if (arg == "--eagm") {
	if (strcmp(argv[i+1],"global") == 0) {
	  eagms.push_back(global);
	} else if (strcmp(argv[i+1],"node") == 0) {
	  eagms.push_back(node);
	} else if (strcmp(argv[i+1],"numa") == 0) {
	  eagms.push_back(numa);
     	} else if (strcmp(argv[i+1],"thread") == 0) {
	  eagms.push_back(thread);
	} else {
	  std::cout << "[ERROR] Invalid EAGM. Available : global, node, numa, thread" << std::endl;
	  return false;
	}
      }

      if (arg == "--sv") {
	// TODO
      }

      if (arg == "--id-distribution") {
	if (strcmp(argv[i+1],"vertical") == 0)
	  id_distribution = vertical;
	else if (strcmp(argv[i+1],"horizontal") == 0)
	  id_distribution = horizontal;
	else {
	  std::cout << "Invalid id distribution type. Available types are vertical and horizontal" 
		    << std::endl;
	  return false;
	}
      }

      if (arg == "--verify") {
	verify = true;
      }
    }
  }

  void print() {
    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "verify : " << verify << std::endl;
  }

  const std::vector<cc_instance_params>&
  get_instance_params() {
    if (instance_params.empty()) {

      if (deltas.empty()) {
	deltas.push_back(10);
      }

      // at least one agm must be specified
      BOOST_FOREACH(int d, deltas) {
	BOOST_FOREACH(cc_agm a, agms) {
	  if (!eagms.empty()) {
	    BOOST_FOREACH(cc_eagm ea, eagms) {
	      instance_params.push_back(cc_instance_params(threads,
                                                           a, 
							   ea,
							   undefined_sv,
							   id_distribution, 
							   verify, 
							   d));
	    }
	  } else {
	    // eagms not specified -- assume chaotic
	    instance_params.push_back(cc_instance_params(threads,
                                                         a, 
							 global, 
							 undefined_sv,
							 id_distribution, 
							 verify, 
							 d));
	  }
	}
      }

      // TODO
      /*BOOST_FOREACH (cc_sv sv, svs) {
	instance_params.push_back(cc_instance_params(threads,
						     undefined_agm, 
						     undefined_eagm, 
						     sv,
						     id_distribution, 
						     verify, 
						     10));

						     }*/

      return instance_params;
    }
  }
};

class CCExecutor {

private:
  template <typename Graph, typename ComponentMap>
  bool verify_cc(Graph& g, ComponentMap& components) {

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
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
    return true;
  }


public:

  CCExecutor() {
  }

  template <typename Graph, 
	    typename Algorithm,
	    typename ComponentsMap,
	    typename graph_create_params>
  time_type
  run_algorithm(amplusplus::transport& trans, 
		amplusplus::transport& barrier_trans, 
		Graph& g, 
		Algorithm& D,
		ComponentsMap& components,
		graph_create_params& gparam,
		instance_params& runtime_param,
		cc_instance_params& algoparams) {
#ifdef CLONE
    amplusplus::transport trans = trans.clone(); // Clone transport for this run
#endif

    trans.set_nthreads(1);
    
    { amplusplus::scoped_epoch epoch(barrier_trans); }

#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS  
    epoch_times.clear();
    clear_buffer_stats();
#endif
    
    // Many threads now
    trans.set_nthreads(runtime_param.threads);

    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_param.threads - 1]);
    for (int i = 0; i < runtime_param.threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < runtime_param.threads - 1; ++i)
      threads[i].join();

    time_type start = D.get_start_time();    
    time_type end = D.get_end_time();

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

    // Verification
    if (algoparams.verify) {
      verify_cc(g, components); 
    }

    time_type exec_time = end-start; 
    algoparams.work_stats.reduce_stats(exec_time); 

    // Stats
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    print_and_clear_epoch_times();
    print_buffer_stats();
#endif


    if (trans.rank() == 0) {
      std::cout << "[INFO] Time : " << (end - start) 
		<< " sec." << std::endl;
    }

    return exec_time;
  }


  template <typename Graph, typename MessageGenerator,
	    typename graph_create_params,
            typename id_dist_t>
  time_type execute_with_id_dist(Graph& g, 
				 amplusplus::transport& trans, 
				 MessageGenerator& msg_gen,
				 graph_create_params& gparams,
				 instance_params& runtime_params,
				 cc_instance_params& _ccparams,
				 id_dist_t & dist) { 

    amplusplus::transport barrier_trans = trans.clone();

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, weight_type WeightedEdge::*>::type WeightMap;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

    // Components map
    std::vector<Vertex> componentS(boost::num_vertices(g), std::numeric_limits<Vertex>::max());
    typedef boost::iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap>  
      ComponentsMap;
    ComponentsMap components(componentS.begin(), get(boost::vertex_index, g));

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      put(components, v, v);
    }

    time_type t = -1;
    trans.set_nthreads(runtime_params.threads);

    if (_ccparams.agm == chaotic) {
      if (_ccparams.eagm == thread) {
	boost::graph::distributed::data_driven_cc<Graph, 
						  ComponentsMap,
						  id_dist_t,
						  agm_work_stats,
						  THREAD_PRIORITY_Q_GEN,
						  MessageGenerator>
	  D(g, 
	    trans,
	    sched_getcpu(),
	    components, 
	    dist,
	    _ccparams.work_stats,
	    msg_gen, 
	    runtime_params.flush);

	t = run_algorithm(trans,
			  barrier_trans,
			  g,
			  D,
			  components,
			  gparams,
			  runtime_params,
			  _ccparams);


      }

      if (_ccparams.eagm == global) {
	boost::graph::distributed::chaotic_cc<Graph, 
					      ComponentsMap,
					      id_dist_t,
					      agm_work_stats,
					      MessageGenerator>
	  D(g, 
	    trans,  
	    sched_getcpu(),
	    components, 
	    dist,
	    _ccparams.work_stats, 
	    msg_gen);

	t = run_algorithm(trans,
			  barrier_trans,
			  g,
			  D,
			  components,
			  gparams,
			  runtime_params,
			  _ccparams);


      }
    } // end of chaotic agm

    if (_ccparams.agm == level) {
      if (_ccparams.eagm == global) {
	boost::graph::distributed::level_sync_cc<Graph, 
						 ComponentsMap,
						 id_dist_t,
						 agm_work_stats,
						 MessageGenerator>
	  D(g, 
	    trans,  
	    sched_getcpu(),
	    components, 
	    dist,
	    _ccparams.work_stats, 
	    msg_gen);

	t = run_algorithm(trans,
			  barrier_trans,
			  g,
			  D,
			  components,
			  gparams,
			  runtime_params,
			  _ccparams);

      }
    }

    if (_ccparams.agm == delta) {
      if (_ccparams.eagm == global) {
	boost::graph::distributed::delta_stepping_cc<Graph, 
						     ComponentsMap,
						     id_dist_t,
						     agm_work_stats,
						     MessageGenerator>
	  D(g, 
	    trans,  
	    sched_getcpu(),
	    components, 
	    _ccparams.delta_val,
	    dist,
	    _ccparams.work_stats, 
	    msg_gen);

	t = run_algorithm(trans,
			  barrier_trans,
			  g,
			  D,
			  components,
			  gparams,
			  runtime_params,
			  _ccparams);

      }

    }

    if (t == -1)
      std::cout << "[ERROR] No algorithm ran ..." << std::endl;

    return t;
  }


public:
  template <typename Graph, typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       cc_instance_params& _ccparams) { 

    amplusplus::transport barrier_trans = trans.clone();

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, weight_type WeightedEdge::*>::type WeightMap;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

    // Components map
    std::vector<Vertex> componentS(boost::num_vertices(g), std::numeric_limits<Vertex>::max());
    typedef boost::iterator_property_map<typename std::vector<Vertex>::iterator, VertexIndexMap>  
      ComponentsMap;
    ComponentsMap components(componentS.begin(), get(boost::vertex_index, g));

    time_type t = -1;
    trans.set_nthreads(runtime_params.threads);


    if (_ccparams.id_distribution == horizontal) {
      typedef row_id_distribution<Graph> id_dist_t;
      id_dist_t dist(g, trans.size());
      return execute_with_id_dist(g,
				  trans,
				  msg_gen,
				  gparams,
				  runtime_params,
				  _ccparams,
				  dist);
    } else if (_ccparams.id_distribution == vertical) {
      typedef block_id_distribution<Graph> id_dist_t;
      id_dist_t dist(g, gparams.n);
      return execute_with_id_dist(g,
				  trans,
				  msg_gen,
				  gparams,
				  runtime_params,
				  _ccparams,
				  dist);

    } else {
      std::cout << "[ERROR] Unsupported distribution type. Only supports horizontal and vertical"
		<< std::endl;
      return -1;
    }
  }

};


int main(int argc, char* argv[]) {
  std::cout << "printing core id for process ..." << std::endl;
  print_core_id();

  executor<CCExecutor, cc_params, cc_instance_params> cc_executor;
  cc_executor.execute(argc, argv);  
  amplusplus::clear_mpi_datatype_map();
}
