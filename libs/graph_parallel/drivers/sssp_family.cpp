// Copyright (C) 2017 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== SSSP (AGM/EAGM) Hand-coded Algortihms================//
// Driver for SSSP family of algorithms.
//===========================================================//

#include "common/sssp_params.hpp"

#define LIB_CDS 1

#include <boost/graph/distributed/thread_pq_def.hpp>
//delta eagms
#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths_node.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths_numa.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths_thread.hpp>
//kla eagms
#include <boost/graph/distributed/kla_sssp_buffer.hpp>
#include <boost/graph/distributed/kla_sssp_numa.hpp>
#include <boost/graph/distributed/kla_sssp_node.hpp>
#include <boost/graph/distributed/kla_sssp_thread.hpp>
//chaotic agms
#include <boost/graph/distributed/distributed_control.hpp>
#include <boost/graph/distributed/distributed_control_node.hpp>
#include <boost/graph/distributed/distributed_control_chaotic.hpp>

#define NODE_PRIORITY_Q_GEN boost::graph::distributed::node_priority_queue_gen
#define NUMA_PRIORITY_Q_GEN boost::graph::distributed::numa_priority_queue_gen
#define THREAD_PRIORITY_Q_GEN boost::graph::distributed::thread_priority_queue_gen

class SSSPExecutor {

private:
  template <typename Graph, typename DistMap, typename WeightMap>
  bool verify_sssp(amplusplus::transport& trans,  Graph& g, DistMap& distance, WeightMap& weight) {
    distance.set_consistency_model(boost::parallel::cm_forward);
    distance.set_max_ghost_cells(0);

    if (trans.rank() == 0) std::cout<<"[INFO] Verifying results......";

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
	if (get(distance, target(e, g)) > 
	    boost::closed_plus<weight_type>()(get(distance, source(e, g)), get(weight, e)))
	  std::cout << get(get(vertex_local, g), source(e, g)) << "@" 
		    << get(get(vertex_owner, g), source(e, g)) << "->"
		    << get(get(vertex_local, g), target(e, g)) << "@" 
		    << get(get(vertex_owner, g), target(e, g)) << "  weight = "
		    << get(weight, e)
		    << "  distance(" << get(get(vertex_local, g), source(e, g)) << "@" 
		    << get(get(vertex_owner, g), source(e, g))
		    << ") = " << get(distance, source(e, g)) << "  distance(" 
		    << get(get(vertex_local, g), target(e, g)) << "@" 
		    << get(get(vertex_owner, g), target(e, g)) << ") = " 
		    << get(distance, target(e, g)) << std::endl;
#else
	if(get(distance, target(e, g)) > boost::closed_plus<weight_type>()(get(distance, v), get(weight, e))) std::abort();
#endif
      }
    }

    if (trans.rank() == 0) std::cout << "Verified." << std::endl;
    distance.clear(); // Clear memory used by ghost cells

    return true;
  }


  template <typename Graph, 
	    typename DistanceMap>
  uint64_t count_visited_edges(Graph& g, DistanceMap& distance) {

    uint64_t vedges = 0;
    weight_type dist;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      dist = boost::get(distance, v);
      if(dist < std::numeric_limits<weight_type>::max()) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  // check whether this is a self loop; if self loop dont count it
	  if (target(e,g) == v)
	    continue;

	  vedges+=1;
	}
      }
    }

    uint64_t totaledges = 0;
    MPI_Allreduce(&vedges, &totaledges, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    // Graph is undirected, Therefore, divide by two
    return (totaledges/2);
  }

  template <typename Graph, 
	    typename DistanceMap>
  uint64_t count_visited_vertices(Graph& g, DistanceMap& distance) {
    uint64_t vvertices = 0;
    weight_type dist;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      dist = boost::get(distance, v);
      if(dist < std::numeric_limits<weight_type>::max()) {
	vvertices+=1;
      }
    }

    uint64_t totalverts = 0;
    MPI_Allreduce(&vvertices, &totalverts, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    return totalverts;
  }

  void init_cds() {
#ifdef LIB_CDS
    // LIBCDS Stuff goes here
    if (_RANK == 0)
      std::cout << "[INFO] Initializing LIBCDS .... " << std::endl;

    cds::Initialize();
    // Initialize Hazard Pointer singleton
    cds::gc::HP hpGC;
    // Attach for the main thread
    cds::threading::Manager::attachThread();
#endif
  }

  void terminate_cds() {
    if (_RANK == 0)
      std::cout << "[INFO] Termminating LIBCDS .... " << std::endl;

    cds::Terminate();
  }

  

public:

  SSSPExecutor() {
  }

  template <typename Graph, 
	    typename Algorithm,
	    typename WeightMap,
	    typename DistanceMap,
	    typename graph_create_params>
  time_type
  run_algorithm(amplusplus::transport& trans, 
		amplusplus::transport& barrier_trans, 
		Graph& g, 
		Algorithm& D,
		WeightMap& weight,
		DistanceMap& distance,
		graph_create_params& gparam,
		instance_params& runtime_param,
		sssp_instance_params& ssspparam) {
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
    
    time_type end = get_time();
    time_type start = D.get_start_time();

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

    // Verification
    if (ssspparam.verify) {
      verify_sssp(trans, g, distance, weight); 
    }

    uint64_t total_vertices = count_visited_vertices(g, distance);
    if (total_vertices < ssspparam.visited_threshold) { 
      ssspparam.work_stats.reset();
      return -1.;
    }

    time_type exec_time = end-start;
    ssspparam.work_stats.reduce_stats(exec_time);

    // Stats
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS
    print_and_clear_epoch_times();
    print_buffer_stats();
#endif

    uint64_t edges_in_comp =  count_visited_edges(g, distance);

    if (trans.rank() == 0) {
      std::cout << "[INFO] Edge traversals should be : " << edges_in_comp
		<< ", vertices visited : " << total_vertices
		<< ", Time : " << exec_time 
		<< " sec." << std::endl;
    }

    return exec_time;
  }


public:
  template <typename Graph, typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       sssp_instance_params& sssp_params) { 

    amplusplus::transport barrier_trans = trans.clone();

    // g must be weighted graph
    typedef boost::compressed_sparse_row_graph<boost::directedS, 
					       boost::no_property, 
					       WeightedEdge, 
					       boost::no_property, 
					       boost::distributedS<unsigned long long> > WGraph;
    //    static_assert(std::is_same<Graph,WGraph>::value, "Graph must be a weighted graph");
    // bad way of casting!
    WGraph& weightedgraph = (WGraph&)(g);

    typedef typename boost::graph_traits<WGraph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<WGraph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<WGraph, weight_type WeightedEdge::*>::type WeightMap;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

    WeightMap weight = boost::get(&WeightedEdge::weight, weightedgraph);

    // Distance map
    std::vector<weight_type> distanceS(boost::num_vertices(weightedgraph), std::numeric_limits<weight_type>::max());
    typedef boost::iterator_property_map<std::vector<weight_type>::iterator, VertexIndexMap>  DistanceMap;
    DistanceMap distance(distanceS.begin(), get(boost::vertex_index, g));

    typedef typename boost::property_traits<WeightMap>::value_type Dist;

    time_type t = -1;
    while(t == -1) {

      Vertex current_source = boost::vertex(sssp_params.get_source(gparams.n), weightedgraph);    
      // we need this as some algorithms retrive threads from the transport
      // we need this here incase if we hit a t=-1
      trans.set_nthreads(runtime_params.threads);

      if (_RANK == 0)
	std::cout << "[INFO] Executing SSSP for source : " << current_source << std::endl;

      if (sssp_params.agm == delta) {
	if (sssp_params.eagm == global) {

	  boost::graph::distributed::delta_stepping_shortest_paths<WGraph, 
								   DistanceMap, 
								   WeightMap,
								   agm_work_stats,
								   append_buffer<Vertex, 
										 10u>, 
								   MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (Dist)sssp_params.delta_val, 
	      sched_getcpu(), 
	      sssp_params.work_stats,
	      msg_gen);

	  D.set_source(current_source);

	  if (sssp_params.level_sync)
	    D.set_level_sync();

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
        }

	if (sssp_params.eagm == node) {

	  init_cds();
	  boost::graph::distributed::delta_stepping_shortest_paths_node<WGraph, 
									DistanceMap, 
									WeightMap,
									agm_work_stats,
									NODE_PRIORITY_Q_GEN, 
									MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (Dist)sssp_params.delta_val, 
	      sched_getcpu(),
	      sssp_params.work_stats,
	      runtime_params.flush, 
	      msg_gen);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
	  terminate_cds();
	} 

	if (sssp_params.eagm == numa) {
	  init_cds();
	  boost::graph::distributed::delta_stepping_shortest_paths_numa<WGraph, 
									DistanceMap, 
									WeightMap,
									agm_work_stats,
									NUMA_PRIORITY_Q_GEN, 
									MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (Dist)sssp_params.delta_val, 
	      sched_getcpu(),
	      sssp_params.work_stats,
	      runtime_params.flush, 
	      msg_gen);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);

	  terminate_cds();
	} 

	if (sssp_params.eagm == thread) {

	  boost::graph::distributed::delta_stepping_shortest_paths_thread<WGraph, 
									  DistanceMap, 
									  WeightMap,
									  agm_work_stats,
									  THREAD_PRIORITY_Q_GEN, 
									  MessageGenerator>
	    D(weightedgraph, distance, weight, trans, (Dist)sssp_params.delta_val, 
	      sched_getcpu(), 
	      sssp_params.work_stats,
	      runtime_params.flush,  msg_gen);
	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
	} 
      } // end of delta EAGMs


      if (sssp_params.agm == kla) {
	if (sssp_params.eagm == global) {

	  boost::graph::distributed::kla_shortest_paths_buffer<WGraph, 
							       DistanceMap, 
							       WeightMap,
							       agm_work_stats,
							       append_buffer<std::pair<Vertex, std::pair<Dist, size_t> >, 10u>,
							       MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (size_t)sssp_params.k_val,
	      sched_getcpu(),
	      sssp_params.work_stats, 
	      msg_gen);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
	} 

	if (sssp_params.eagm == node) {
	  init_cds();
	  boost::graph::distributed::kla_shortest_paths_node<WGraph, 
							     DistanceMap, 
							     WeightMap,
							     agm_work_stats,
							     NODE_PRIORITY_Q_GEN,
							     MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (size_t)sssp_params.k_val,
	      sched_getcpu(),
	      sssp_params.work_stats, 
              runtime_params.flush,
	      msg_gen);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);

          terminate_cds();
        }

	if (sssp_params.eagm == numa) {
	  init_cds();
	  boost::graph::distributed::kla_shortest_paths_numa<WGraph, 
							     DistanceMap, 
							     WeightMap,
							     agm_work_stats,
							     NUMA_PRIORITY_Q_GEN,
							     MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (size_t)sssp_params.k_val,
	      sched_getcpu(),
	      sssp_params.work_stats, 
              runtime_params.flush,
	      msg_gen);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
          terminate_cds();
        }

	if (sssp_params.eagm == thread) {

	  boost::graph::distributed::kla_shortest_paths_thread<WGraph, 
							       DistanceMap, 
							       WeightMap,
							       agm_work_stats,
							       THREAD_PRIORITY_Q_GEN,
							       MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      (size_t)sssp_params.k_val,
	      sched_getcpu(),
	      sssp_params.work_stats, 
              runtime_params.flush,
	      msg_gen);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
        }
      } // end of kla EAGMs



      if (sssp_params.agm == chaotic) {
	if (sssp_params.eagm == global) {

	  boost::graph::distributed::distributed_control_chaotic<WGraph, 
								 DistanceMap, 
								 WeightMap,
								 agm_work_stats,
								 MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      sched_getcpu(),
	      sssp_params.work_stats, 
	      msg_gen,
              runtime_params.flush);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
	} 

	if (sssp_params.eagm == node) {
	  init_cds();
	  boost::graph::distributed::distributed_control_node<WGraph, 
							      DistanceMap, 
							      WeightMap,
							      agm_work_stats,
							      NODE_PRIORITY_Q_GEN,
							      MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      sched_getcpu(),
	      sssp_params.work_stats, 
	      msg_gen,
              runtime_params.flush,
              false);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
          terminate_cds();
        }

	if (sssp_params.eagm == numa) {
	  init_cds();
	  boost::graph::distributed::distributed_control_node<WGraph, 
							      DistanceMap, 
							      WeightMap,
							      agm_work_stats,
							      NUMA_PRIORITY_Q_GEN,
							      MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      sched_getcpu(),
	      sssp_params.work_stats, 
	      msg_gen,
              runtime_params.flush,
              true); //NUMA

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
          terminate_cds();

        }

	if (sssp_params.eagm == thread) {

	  boost::graph::distributed::distributed_control<WGraph, 
							 DistanceMap, 
							 WeightMap,
							 agm_work_stats,
							 boost::graph::distributed::default_priority_queue_gen,
							 MessageGenerator>
	    D(weightedgraph, 
	      distance, 
	      weight, 
	      trans, 
	      sched_getcpu(),
	      sssp_params.work_stats, 
              msg_gen,
              runtime_params.flush);

	  D.set_source(current_source);

	  t = run_algorithm(trans,
			    barrier_trans,
			    weightedgraph,
			    D,
			    weight,
			    distance,
			    gparams,
			    runtime_params,
			    sssp_params);
	}
      } // end of Chaotic EAGMs


      if (t == -1) {
	if (_RANK == 0)
	  std::cout << "[WARNING] Execution did not traversed threshold amount of vertices."
		    << " Threshold : " << sssp_params.visited_threshold 
		    << ". Continuing..."
		    << std::endl;
      } else
	return t;
    } // end of while(t == -1) 

  }

};


int main(int argc, char* argv[]) {
  std::cout << "printing core id for process ..." << std::endl;
  print_core_id();

  executor<SSSPExecutor, sssp_params, sssp_instance_params> sssp_executor;
  sssp_executor.execute(argc, argv);  
  amplusplus::clear_mpi_datatype_map();
}
