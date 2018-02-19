#ifndef __AGM_SSSP_HPP
#define __AGM_SSSP_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/model/agm.hpp>
#include <boost/graph/agm/runtime/runtime.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>

#define VISITED_THRESHOLD 100

namespace boost { namespace graph { namespace agm {
      
// definition of the SSSP work item set
typedef uint32_t Distance;
typedef uint32_t Weight;

template<typename Graph>
class sssp_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;

  // work item set definition
  // Every work item must have the time
  // destination, source, distance
  typedef std::tuple<Vertex, Distance> WorkItem;

  // For level ordering, redifining the work item
  typedef int Level;
  // work item set definition
  // Every work item must have the time
  // destination, distance
  typedef std::tuple<Vertex, Distance, Level> LevelWorkItem;

  //===================================================================================
  // Pair of processing functions : post order pf and pre-order pf
  //===================================================================================
  // the processing function  
  template<typename WeightState,
           typename DistanceState>
  struct post_order_sssp_pf {
    
  public:
    post_order_sssp_pf(const Graph& _rg,
                       WeightState& _vweight,
                       DistanceState& _dist,
                       agm_work_stats& _sr) : g(_rg),
                                              vdistance(_dist),
                                              vweight(_vweight),
                                              stats(_sr){}
    
  private:
    const Graph& g;
    DistanceState& vdistance;
    WeightState& vweight;
    agm_work_stats& stats;    

  public:
    template<typename outset>
    void operator()(const WorkItem& w,
                    int tid,
                    outset& out) {
      Vertex v = std::get<0>(w);
      Distance d = std::get<1>(w);

      if (vdistance[v] == d) {
        BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
          Vertex u = boost::target(e, g);
          Distance we = boost::get(vweight, e);
#ifdef PBGL2_PRINT_WORK_STATS
          stats.increment_edges(tid);
#endif          
          WorkItem generated(u, (d+we));
          out.push(generated, tid);
        }
      }
#ifdef PBGL2_PRINT_WORK_STATS
      else {
        stats.increment_invalidated_cancels(tid);
      }
#endif                
    }

    // for level work items
    template<typename outset>
    void operator()(const LevelWorkItem& w,
                    int tid,
                    outset& out) {
      Vertex v = std::get<0>(w);
      Distance d = std::get<1>(w);
      Level level = std::get<2>(w);            

      if (vdistance[v] == d) {
        BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
          Vertex u = boost::target(e, g);
          Distance we = boost::get(vweight, e);
#ifdef PBGL2_PRINT_WORK_STATS
          stats.increment_edges(tid);
#endif          
          LevelWorkItem generated(u, (d+we), (level+1));
          out.push(generated, tid);
        }
      }
#ifdef PBGL2_PRINT_WORK_STATS
      else {
        stats.increment_invalidated_cancels(tid);
      }
#endif      
    }    
  };  
  
  // the pre-order processing function  
  template<typename WeightState,
           typename DistanceState>
  struct pre_order_sssp_pf {

  public:
    pre_order_sssp_pf(const Graph& _rg,
            WeightState& _weight,
            DistanceState& _dist,
            agm_work_stats& _sr) : g(_rg),
                                   vweight(_weight),
                                   vdistance(_dist),
                                   stats(_sr){}

  private:
    const Graph& g;
    WeightState& vweight;        
    DistanceState& vdistance;
    agm_work_stats& stats;

  public:

    template<typename work_item, typename buckets>
    void operator()(const work_item& wi,
                    int tid,
		    buckets& outset) {
      Vertex v = std::get<0>(wi);
      Distance d = std::get<1>(wi);

      Distance old_dist = vdistance[v], last_old_dist;

      while (d < old_dist) {
        last_old_dist = old_dist;
        old_dist = boost::parallel::val_compare_and_swap(&vdistance[v], old_dist, d);
        if (last_old_dist == old_dist) {
#ifdef PBGL2_PRINT_WORK_STATS        
          if (old_dist < std::numeric_limits<weight_type>::max()) {
            stats.increment_invalidated(tid);
          } else
            stats.increment_useful(tid);            
#endif
          outset.push(wi, tid);          
          return;
        }
      }
      
#ifdef PBGL2_PRINT_WORK_STATS                
      stats.increment_rejected(tid);
#endif      
    }
  };

  //=====================================================================
  // The general processing function that can be used as the
  // pre-order processing function or post-order processing
  // function.
  //=====================================================================
  // the processing function  
  template<typename WeightState,
           typename DistanceState>
  struct sssp_pf {

  public:
    sssp_pf(const Graph& _rg,
            WeightState& _weight,
            DistanceState& _dist,
            agm_work_stats& _sr) : g(_rg),
                                   vweight(_weight),
                                   vdistance(_dist),
                                   stats(_sr){}

  private:
    const Graph& g;
    WeightState& vweight;        
    DistanceState& vdistance;
    agm_work_stats& stats;

  public:

    template<typename buckets>
    void operator()(const WorkItem& wi,
                    int tid,
		    buckets& outset) {

      Vertex v = std::get<0>(wi);
      Distance d = std::get<1>(wi);

      Distance old_dist = vdistance[v], last_old_dist;

      while (d < old_dist) {
        last_old_dist = old_dist;
        old_dist = boost::parallel::val_compare_and_swap(&vdistance[v], old_dist, d);
        if (last_old_dist == old_dist) {
#ifdef PBGL2_PRINT_WORK_STATS        
          if (old_dist < std::numeric_limits<weight_type>::max()) {
            stats.increment_invalidated(tid);
          } else
            stats.increment_useful(tid);            
#endif
          
          BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
            Vertex u = boost::target(e, g);
            Distance we = boost::get(vweight, e);
#ifdef PBGL2_PRINT_WORK_STATS
            stats.increment_edges(tid);
#endif
            WorkItem generated(u, (d+we));
            outset.push(generated, tid);
          }

          return;
        }
      }
      
#ifdef PBGL2_PRINT_WORK_STATS                
      stats.increment_rejected(tid);
#endif      
    }


    // For level ordering
    template<typename buckets>
    void operator()(const LevelWorkItem& wi,
                    int tid,
		    buckets& outset) {
      
      Vertex v = std::get<0>(wi);
      Distance d = std::get<1>(wi);
      Level l = std::get<2>(wi);      

      Distance old_dist = vdistance[v], last_old_dist;
      
      while (d < old_dist) {
        last_old_dist = old_dist;
        old_dist = boost::parallel::val_compare_and_swap(&vdistance[v], old_dist, d);
        if (last_old_dist == old_dist) {
#ifdef PBGL2_PRINT_WORK_STATS        
          if (old_dist < std::numeric_limits<weight_type>::max()) {
            stats.increment_invalidated(tid);
          } else
            stats.increment_useful(tid);            
#endif
          
          BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
            Vertex u = boost::target(e, g);
            Distance we = boost::get(vweight, e);
#ifdef PBGL2_PRINT_WORK_STATS
            stats.increment_edges(tid);
#endif            
            LevelWorkItem generated(u, (d+we), (l+1));
            outset.push(generated, tid);
          }

          return;
        }
      }
      
#ifdef PBGL2_PRINT_WORK_STATS                
      stats.increment_rejected(tid);
#endif      
    }
  };
  

public:

  template <typename DistMap, typename WeightMap>
  bool verify(const Graph& g, DistMap& distance, WeightMap& weight) {
    distance.set_consistency_model(boost::parallel::cm_forward);
    distance.set_max_ghost_cells(0);

    { amplusplus::scoped_epoch epoch(g.transport()); }
    
    if (g.transport().rank() == 0) std::cout<<"[INFO] Verifying results......";

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
        //#ifdef 1 //PRINT_DEBUG
	if (get(distance, target(e, g)) > 
	    boost::closed_plus<Distance>()(get(distance, v), get(weight, e))) {
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
          
          std::cout << "[ERROR] Verification failed." << std::endl;
          assert(false);
        }
      }
    }

    if (g.transport().rank() == 0) std::cout << "Verified." << std::endl;
    distance.clear(); // Clear memory used by ghost cells

    return true;
  }

  template <typename DistanceMap>
  uint64_t count_visited_vertices(Graph& g, DistanceMap& distance) {
    uint64_t vvertices = 0;
    weight_type dist;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      dist = boost::get(distance, v);
      if(dist < std::numeric_limits<Distance>::max()) {
	vvertices+=1;
      }
    }

    uint64_t totalverts = 0;
    MPI_Allreduce(&vvertices, &totalverts, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    return totalverts;
  }
  

  template<typename WeightMap,
           typename DistMap,
           typename RuntimeModelGen,
           typename EAGMConfig,
           typename WorkItemType>
  time_type execute_split_pf(Graph& g,
                             WeightMap& weights,
                             DistMap& distance_state,
                             RuntimeModelGen rtmodelgen,
                             EAGMConfig& config,
                             WorkItemType sw,
                             instance_params& runtime_params,
                             agm_work_stats& sr, 
                             bool _verify=true) {
    
    info("Setting the initial work item set ...");
    // Initial work item set
    typedef append_buffer<WorkItemType, 10u> InitialWorkItems;
    InitialWorkItems initial;
    auto s = std::get<0>(sw);
    if (get(get(vertex_owner, g), s) == _RANK) {      
      distance_state[s] = 0;
      initial.push_back(sw);
      
#ifdef PBGL2_PRINT_WORK_STATS        
      sr.increment_useful(0);            
#endif
      
    }
    
    // Processing funcion
    typedef pre_order_sssp_pf<WeightMap, DistMap> ProcessingFunction;
    ProcessingFunction pf(g, weights, distance_state, sr);

    //Send predicate
    typedef post_order_sssp_pf<WeightMap, DistMap> PostOrderProcessingFunction;
    PostOrderProcessingFunction sendpf(g, weights, distance_state, sr);
    
    // SSSP algorithm
    typedef eagm<Graph,
                 WorkItemType,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen,
                 PostOrderProcessingFunction> sssp_eagm_t;

    sssp_eagm_t ssspalgo(rtmodelgen,
                         config,
                         pf,
                         sendpf,
                         initial);
    
    info("Invoking SSSP algorithm with split processing functions ...");
    time_type elapsed = ssspalgo(runtime_params);

#ifdef PBGL2_PRINT_WORK_STATS
    ssspalgo.print_stats();
#endif    
    
    return elapsed;    
  }
  
  template<typename WeightMap,
           typename DistMap,
           typename RuntimeModelGen,
           typename EAGMConfig,
           typename WorkItemType>
  time_type execute_pre_order_pf(Graph& g,
                                 WeightMap& weights,
                                 DistMap& distance_state,
                                 RuntimeModelGen rtmodelgen,
                                 EAGMConfig& config,
                                 WorkItemType sw,
                                 instance_params& runtime_params,
                                 agm_work_stats& sr, 
                                 bool _verify=true) {
    
    info("Setting the initial work item set ...");
    // Initial work item set
    typedef append_buffer<WorkItemType, 10u> InitialWorkItems;
    InitialWorkItems initial;
    auto s = std::get<0>(sw);
    if (get(get(vertex_owner, g), s) == _RANK) {
      initial.push_back(sw);
    }
    
    // Processing funcion
    typedef sssp_pf<WeightMap, DistMap> ProcessingFunction;
    ProcessingFunction pf(g, weights, distance_state, sr);    
    
    // SSSP algorithm
    typedef eagm<Graph,
                 WorkItemType,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen,
                 EMPTY_PF> sssp_eagm_t;

    sssp_eagm_t ssspalgo(rtmodelgen,
                         config,
                         pf,
                         initial);
    
    info("Invoking SSSP algorithm with preorder processing functions ...");
    time_type elapsed = ssspalgo(runtime_params);

#ifdef PBGL2_PRINT_WORK_STATS
    ssspalgo.print_stats();
#endif    
    
    return elapsed;    
  }  

  
  template<typename WeightMap,
           typename DistMap,
           typename RuntimeModelGen,
           typename EAGMConfig,
           typename WorkItemType>  
  time_type execute_post_order_pf(Graph& g,
                                 WeightMap& weights,
                                 DistMap& distance_state,
                                 RuntimeModelGen rtmodelgen,
                                 EAGMConfig& config,
                                 WorkItemType sw,
                                 instance_params& runtime_params,
                                 agm_work_stats& sr, 
                                 bool _verify=true) {
    
    info("Setting the initial work item set ...");
    // Initial work item set
    typedef append_buffer<WorkItemType, 10u> InitialWorkItems;
    InitialWorkItems initial;
    auto s = std::get<0>(sw);
    if (get(get(vertex_owner, g), s) == _RANK) {
      initial.push_back(sw);
    }
    
    // Processing funcion
    typedef sssp_pf<WeightMap, DistMap> ProcessingFunction;
    ProcessingFunction pf(g, weights, distance_state, sr);
    
    // SSSP algorithm
    typedef eagm<Graph,
                 WorkItemType,
                 EMPTY_PF,
                 EAGMConfig,
                 RuntimeModelGen,
                 ProcessingFunction> sssp_eagm_t;

    sssp_eagm_t ssspalgo(rtmodelgen,
                         config,
                         pf,
                         initial);
    
    info("Invoking SSSP algorithm with post-processing functions ...");
    time_type elapsed = ssspalgo(runtime_params);
    
#ifdef PBGL2_PRINT_WORK_STATS
    ssspalgo.print_stats();
#endif    

    return elapsed;        
  }
  
  template<typename WeightMap,
           typename RuntimeModelGen,
           typename EAGMConfig,
           typename WorkItemType,
           typename agm_param_type>
  time_type execute_eagm_common(Graph& g,
                                WeightMap& weights,
                                RuntimeModelGen rtmodelgen,
                                EAGMConfig& config,
                                WorkItemType sw,
                                agm_param_type& agm_params,
                                instance_params& runtime_params,
                                agm_work_stats& sr, 
                                bool _verify=true) {

    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<Distance> distmap(num_vertices(g), std::numeric_limits<Distance>::max());
    typedef boost::iterator_property_map<typename std::vector<Distance>::iterator, VertexIndexMap> DistMap;
    DistMap distance_state(distmap.begin(), get(boost::vertex_index, g));
    time_type elapsed = -1;

    if (agm_params.pf_mode == agm_pf_splitted) {
      elapsed = execute_split_pf(g,
                                 weights,
                                 distance_state,
                                 rtmodelgen,
                                 config,
                                 sw,
                                 runtime_params,
                                 sr,
                                 _verify);
    } else if (agm_params.pf_mode == agm_pf_preorder) {
      elapsed = execute_pre_order_pf(g,
                                     weights,
                                     distance_state,
                                     rtmodelgen,
                                     config,
                                     sw,
                                     runtime_params,
                                     sr,
                                     _verify);      
    } else if (agm_params.pf_mode == agm_pf_postorder) {
      elapsed = execute_post_order_pf(g,
                                     weights,
                                     distance_state,
                                     rtmodelgen,
                                     config,
                                     sw,
                                     runtime_params,
                                     sr,
                                     _verify);      
      
    } else {
      error("Invalid processing function invocation mode!");
      assert(false);
    }
    
    if (_verify) {
      verify(g, distance_state, weights);
    }

    //visited vertices threshold
    // if the number if visited vertices are less than
    // visited threshold repeat the run
    auto visited = count_visited_vertices(g, distance_state);
    if (boost::num_vertices(g) > VISITED_THRESHOLD) {
      if (visited < VISITED_THRESHOLD)
        elapsed = -1;
    }

    return elapsed;
  }

  template<typename WeightMap,
           typename RuntimeModelGen,
           typename EAGMConfig,
           typename agm_param_type>
  time_type execute_eagm(Graph& g,
                         WeightMap& weights,
                         RuntimeModelGen rtmodelgen,
                         EAGMConfig& config,
                         agm_param_type& agm_param,
                         Vertex source,
                         instance_params& runtime_params,
                         agm_work_stats& sr, 
                         bool _verify=true) {

    WorkItem w(source, 0);

    return execute_eagm_common(g,
                               weights,
                               rtmodelgen,
                               config,
                               w,
                               agm_param,
                               runtime_params,
                               sr,
                               _verify);
  }  

  template<typename WeightMap,
           typename RuntimeModelGen,
           typename EAGMConfig,
           typename agm_param_type>
  time_type execute_klevel_eagm(Graph& g,
                                WeightMap& weights,
                                RuntimeModelGen rtmodelgen,
                                EAGMConfig& config,
                                agm_param_type& agm_param,
                                Vertex source,
                                instance_params& runtime_params,
                                agm_work_stats& sr, 
                                bool _verify=true) {
    
    LevelWorkItem w(source, 0, 0);

    return execute_eagm_common(g,
                               weights,
                               rtmodelgen,
                               config,
                               w,
                               agm_param,
                               runtime_params,
                               sr,
                               _verify);
  }  

  
};

}}}
#endif
