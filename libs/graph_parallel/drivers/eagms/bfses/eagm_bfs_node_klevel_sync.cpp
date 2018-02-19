#include "../../common/synthetic_generator.hpp"
#include "../../common/graph_loader.hpp"
#include "../../common/parser.hpp"
#include "../../common/executor.hpp"
//#include "common/vertex_permutations.hpp"

#include <boost/graph/util/work_stats.hpp>
#include <boost/graph/agm/model/eagm_buckets.hpp>
#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/algorithms/bfs.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/runtime/ampp_runtime.hpp>

#include <limits.h>

// Static initializations
using general_orderings = boost::graph::agm::ordering_type_map;
general_orderings::ordering_map_t general_orderings::ordering_map = general_orderings::create_ordering_map();

const uint64_t undefined_source = std::numeric_limits<uint64_t>::max();

typedef boost::graph::agm::base_ordering base_ordering_t;

class agm_instance_params {
public:
  int threads;
  int k_val;
  uint64_t source;
  bool verify;
  agm_work_stats work_stats;
  runtime_stats rt_stats;

  
  agm_instance_params(int _threads,
                      int _k,
                      uint64_t s,
                      bool _v) : threads(_threads),
                                 k_val(_k),
                                 source(s),
                                 verify(_v),
                                 work_stats(threads),
                                 rt_stats(threads){}

  void print() {
    std::cout << "source : " << source << std::endl;
    std::cout << "k-val : " << k_val << std::endl;    
    std::cout << "verify : " << verify << std::endl;    
  }


  void print_summary() {
    work_stats.print_summary();    
  }

  std::string get_algorithm() {
    return "AGM-Level-Synchronous";
  }
};


class agm_params {

private:
  std::vector<agm_instance_params> params;
  int threads;
  std::vector<uint64_t> sources;
  std::vector<int> ks;  
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

      if (arg == "--source") {
	sources = extract_params<uint64_t>( argv[i+1] );
      }      

      if (arg == "--k-val") {
	ks = extract_params<int>( argv[i+1] );
      }      
      
      if (arg == "--verify") {
        verify = true;
      }      
    }
    
    return true;
  }

  void print() {
  }

  const std::vector<agm_instance_params>&
  get_instance_params() {
    if (params.size() == 0) {

      if (sources.empty())
	sources.push_back(undefined_source);
      
      if (ks.empty())
	ks.push_back(2);
  
      if (threads == -1) {
	std::cout << "[ERROR] Number of threads must be specified using --threads ..." << std::endl;
	assert(false);
      }

      BOOST_FOREACH (uint64_t s, sources) {
        BOOST_FOREACH (int k, ks) {
          agm_instance_params ainst(threads, k, s, verify);
          params.push_back(ainst);
        }
      }

    }
    
    return params;
  }
};

class AGMBFSNodeKLevelSyncExecutor {

public:
  template <typename Graph, typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(const Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       agm_instance_params& agm_params) {


    boost::graph::agm::bfs_family<Graph> bfses;
    CHAOTIC_ORDERING_T ch;
    KLEVEL_ORDERING_T klevel(agm_params.k_val);    
    auto config = boost::graph::agm::create_eagm_config(ch,
                                                        klevel,
                                                        ch,
                                                        ch);

    info("Running with following EAGM configuration;");
    config.print();

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

    time_type exec_time = bfses.execute_eagm(g,
                                             rtgen,
                                             config,
                                             runtime_params,
                                             agm_params.work_stats,
                                             agm_params.verify);
    
    agm_params.work_stats.reduce_stats(exec_time);
    return exec_time;
    
  }
};

int main(int argc, char* argv[]) {

  boost::graph::agm::ordering_type_map::init();
  info("Starting AGM ...");
  std::cout << "[INFO] Printing core id for process ..." << std::endl;
  print_core_id();

  executor<AGMBFSNodeKLevelSyncExecutor, agm_params, agm_instance_params> agm_executor;
  agm_executor.execute(argc, argv);
  boost::graph::agm::ordering_type_map::clear();   
}
