#include "common/sssp_params.hpp"

#include <boost/graph/agm/util/bucket.hpp>
#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/algorithms/sssp.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>

#include <limits.h>

class AGMSSSPExecutor {

private:

  template <typename Graph,
            typename MessageGenerator,
            typename WeightMap,
            typename OwnerMap,
            typename SSSPFamily,
	    typename graph_create_params>
  time_type execute_agm(const Graph& g, 
                        amplusplus::transport& trans, 
                        MessageGenerator& msg_gen,
                        WeightMap& weight,
                        OwnerMap& owner,
                        SSSPFamily& ssspes,
                        graph_create_params& gparams,
                        instance_params& runtime_params,
                        sssp_instance_params& agm_params) {
    // runtime
    boost::graph::agm::runtime_wrapper_gen<MessageGenerator,
                                           OwnerMap> rtgen(trans,
                                                           msg_gen,
                                                           runtime_params,
                                                           owner,
                                                           sched_getcpu());
    time_type exec_time = -1;
    
    if (agm_params.agm == sssp_agm::chaotic) {
      exec_time = ssspes.execute(g,
                                 weight,
                                 rtgen,
                                 boost::graph::agm::chaotic(),
                                 runtime_params,
                                 agm_params.work_stats,
                                 agm_params.verify);
      
    } else if (agm_params.agm == sssp_agm::dijkstra) {
      exec_time = ssspes.execute(g,
                                 weight,
                                 rtgen,
                                 boost::graph::agm::dijkstra<1>(),
                                 runtime_params,
                                 agm_params.work_stats,
                                 agm_params.verify);
      
    } else if (agm_params.agm == sssp_agm::delta) {
      exec_time = ssspes.execute(g,
                                 weight,
                                 rtgen,
                                 boost::graph::agm::delta_ord<1, weight_type>(agm_params.delta_val),
                                 runtime_params,
                                 agm_params.work_stats,
                                 agm_params.verify);
    } else {
      error("Unknown AGM algorithm : ", agm_params.agm);
      assert(false);
    }

    agm_params.work_stats.reduce_stats(exec_time);
    return exec_time;
  }

  template<typename T1,
           typename T2,
           typename T3,
           typename T4>
  boost::graph::agm::eagm_configs<T1,
                                  T2,
                                  T3,
                                  T4>
  buildEAGMConfig(T1 _t1,
                  T2 _t2,
                  T3 _t3,
                  T4 _t4,
                  instance_params& rt_params,
                  int ranks){
    
    typedef boost::graph::agm::eagm_configs<T1,
                                            T2,
                                            T3,
                                            T4> eagm_config_t;
    eagm_config_t config(_t1,
                         _t2,
                         _t3,
                         _t4);
    config.ranks = ranks;
    config.threads = rt_params.threads;
    config.optimize();

    return config;
  }
  
  template <typename Graph,
            typename MessageGenerator,
            typename WeightMap,
            typename OwnerMap,
            typename SSSPFamily,
	    typename graph_create_params>
  time_type execute_eagm(const Graph& g, 
                        amplusplus::transport& trans, 
                        MessageGenerator& msg_gen,
                        WeightMap& weight,
                        OwnerMap& owner,
                        SSSPFamily& ssspes,
                        graph_create_params& gparams,
                        instance_params& runtime_params,
                        sssp_instance_params& agm_params) {

    info("inside EAGM execute ...");
    
    // build the eagm config, right now lets hard code DC
    typedef boost::graph::agm::eagm_configs<CHAOTIC_ORDERING_T,
                                            CHAOTIC_ORDERING_T,
                                            CHAOTIC_ORDERING_T,
                                            DIJKSTRA_ORDERING_T> eagm_config_t;

    CHAOTIC_ORDERING_T ch;
    DIJKSTRA_ORDERING_T dj;
    eagm_config_t eagmconfig = buildEAGMConfig(ch, ch, ch, dj, runtime_params, trans.size());
    eagmconfig.print();
    // runtime
    boost::graph::agm::eagm_runtime_wrapper_gen<MessageGenerator,
                                                OwnerMap> rtgen(trans,
                                                                msg_gen,
                                                                runtime_params,
                                                                owner,
                                                                sched_getcpu());
    time_type exec_time = -1;
    
    if (agm_params.agm == sssp_agm::chaotic) {
      exec_time = ssspes.execute_eagm(g,
                                      weight,
                                      rtgen,
                                      eagmconfig,
                                      runtime_params,
                                      agm_params.work_stats,
                                      agm_params.verify);
      
    } else if (agm_params.agm == sssp_agm::dijkstra) {
      exec_time = ssspes.execute_eagm(g,
                                      weight,
                                      rtgen,
                                      eagmconfig,
                                      runtime_params,
                                      agm_params.work_stats,
                                      agm_params.verify);
      
    } else if (agm_params.agm == sssp_agm::delta) {
      exec_time = ssspes.execute_eagm(g,
                                      weight,
                                      rtgen,
                                      eagmconfig,
                                      runtime_params,
                                      agm_params.work_stats,
                                      agm_params.verify);
    } else {
      error("Unknown AGM algorithm : ", agm_params.agm);
      assert(false);
    }

    agm_params.work_stats.reduce_stats(exec_time);
    return exec_time;
  }

  
public:
  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(const Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       sssp_instance_params& agm_params) {
    info("Inside AGM executor ...");

    // g must be weighted graph
    typedef boost::compressed_sparse_row_graph<boost::directedS, 
					       boost::no_property, 
					       WeightedEdge, 
					       boost::no_property, 
					       boost::distributedS<unsigned long long> > WGraph;
    //    static_assert(std::is_same<Graph,WGraph>::value, "Graph must be a weighted graph");
    // bad way of casting!
    WGraph& weightedgraph = (WGraph&)(g);

    typedef typename boost::property_map<WGraph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = boost::get(&WeightedEdge::weight, weightedgraph);


    boost::graph::agm::sssp_family<WGraph> ssspes;
    // distribution
    typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
    OwnerMap owner(boost::get(boost::vertex_owner, g));

    if (agm_params.eagm == global) {
      return execute_agm(weightedgraph, 
                         trans, 
                         msg_gen,
                         weight,
                         owner,
                         ssspes,
                         gparams,
                         runtime_params,
                         agm_params);
    } else {
      return execute_eagm(weightedgraph, 
                          trans, 
                          msg_gen,
                          weight,
                          owner,
                          ssspes,
                          gparams,
                          runtime_params,
                          agm_params);

    }

  
  }

};

int main(int argc, char* argv[]) {

  info("Starting AGM ...");
  std::cout << "[INFO] Printing core id for process ..." << std::endl;
  print_core_id();

  executor<AGMSSSPExecutor, sssp_params, sssp_instance_params> agm_executor;
  agm_executor.execute(argc, argv);  

 
}
