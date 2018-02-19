#ifndef __EAGM_CC_BASE__
#define __EAGM_CC_BASE__
#include "../../common/synthetic_generator.hpp"
#include "../../common/graph_loader.hpp"
#include "../../common/parser.hpp"
#include "../../common/executor.hpp"

#include <boost/graph/util/vertex_permutations.hpp>
#include <boost/graph/util/work_stats.hpp>

#include <boost/graph/agm/model/eagm_buckets.hpp>
#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/algorithms/cc.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/runtime/ampp_runtime.hpp>

#include <limits.h>
typedef boost::graph::agm::base_ordering base_ordering_t;

class agm_instance_params {
public:
  int threads;
  id_distribution_t id_distribution;
  uint64_t delta;
  bool verify;
  agm_work_stats work_stats;
  runtime_stats rt_stats;

  
  agm_instance_params(int _threads,
                      id_distribution_t& _idd,
                      uint64_t _delta,
                      bool _v) : threads(_threads),
                                 id_distribution(_idd),
                                 delta(_delta),
                                 verify(_v),
                                 work_stats(threads),
                                 rt_stats(threads){}

  void print() {
    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;

    std::cout << "delta : " << delta << std::endl;
    std::cout << "verify : " << verify << std::endl;    
  }


  void print_summary() {
    work_stats.print_summary();    
  }

  std::string get_algorithm() {
    return "AGM-CC-Chaotic";
  }
};


class agm_params {

private:
  std::vector<agm_instance_params> params;
  int threads;
  uint64_t delta;
  id_distribution_t id_distribution = horizontal; // default  
  bool verify;
  
public:
  agm_params() : threads(1),
                 verify(false){}

  bool parse(int argc, char* argv[]){
        
    for (int i=0; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--threads") {
	threads = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--verify") {
        verify = true;
      }

      if (arg == "--delta") {
        delta = boost::lexical_cast<uint64_t>( argv[i+1]);
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
      
    }
    
    return true;
  }

  void print() {
  }

  const std::vector<agm_instance_params>&
  get_instance_params() {
    if (params.size() == 0) {
      if (threads == -1) {
	std::cout << "[ERROR] Number of threads must be specified using --threads ..." << std::endl;
	assert(false);
      }

      agm_instance_params ainst(threads, id_distribution, delta, verify);
      params.push_back(ainst);
    }
    
    return params;
  }
};

class AGMCCBaseExecutor {
protected:
  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename IdDistribution,
            typename EAGMConfig>
  time_type execute_with_eagm_config(Graph& g, 
                                     amplusplus::transport& trans, 
                                     MessageGenerator& msg_gen,
                                     graph_create_params& gparams,
                                     instance_params& runtime_params,
                                     agm_instance_params& agm_params,
                                     IdDistribution& iddist,
                                     EAGMConfig& config) {
    
    info("Running with following EAGM configuration;");
    config.print();

    boost::graph::agm::cc_family<Graph, IdDistribution> cces;    
    
    // distribution
    typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
    OwnerMap owner(boost::get(boost::vertex_owner, g));
    
    // runtime
    boost::graph::agm::ampp_runtime_gen<MessageGenerator,
                                        OwnerMap> rtgen(trans,
                                                        msg_gen,
                                                        owner,
                                                        agm_params.rt_stats,
                                                        sched_getcpu(),
                                                        runtime_params);

    time_type exec_time = cces.execute_eagm(g,
                                            rtgen,
                                            config,
                                            runtime_params,
                                            agm_params.work_stats,
                                            iddist,
                                            agm_params.verify);
    
    agm_params.work_stats.reduce_stats(exec_time);
    return exec_time;    
  }

  virtual ~AGMCCBaseExecutor(){}
  
  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename EAGMConfig>
  time_type select_id_distribution(Graph& g, 
                                   amplusplus::transport& trans, 
                                   MessageGenerator& msg_gen,
                                   graph_create_params& gparams,
                                   instance_params& runtime_params,
                                   agm_instance_params& agm_params,
                                   EAGMConfig& config) {

    if (agm_params.id_distribution == horizontal) {
      typedef row_id_distribution<Graph> id_dist_t;
      id_dist_t dist(g, trans.size());
      return execute_with_eagm_config(g,
                                      trans,
                                      msg_gen,
                                      gparams,
                                      runtime_params,
                                      agm_params,
                                      dist,
                                      config);
    } else if (agm_params.id_distribution == vertical) {
      typedef block_id_distribution<Graph> id_dist_t;
      id_dist_t dist(g, gparams.n);
      return execute_with_eagm_config(g,
                                      trans,
                                      msg_gen,
                                      gparams,
                                      runtime_params,
                                      agm_params,
                                      dist,
                                      config);

    } else {
      std::cout << "[ERROR] Unsupported distribution type. Only supports horizontal and vertical"
		<< std::endl;
      return -1;
    }    
  }
};

#endif
