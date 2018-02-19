#include "../../common/synthetic_generator.hpp"
#include "../../common/graph_loader.hpp"
#include "../../common/parser.hpp"
#include "../../common/executor.hpp"
//#include "common/vertex_permutations.hpp"

#include <boost/graph/util/work_stats.hpp>
#include <boost/graph/agm/model/eagm_buckets.hpp>
#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/algorithms/sssp.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/runtime/ampp_runtime.hpp>

#include <limits.h>

// Static initializations
using general_orderings = boost::graph::agm::ordering_type_map;
general_orderings::ordering_map_t general_orderings::ordering_map = general_orderings::create_ordering_map();

typedef boost::graph::agm::base_ordering base_ordering_t;

const uint64_t undefined_source = std::numeric_limits<uint64_t>::max();

class agm_instance_params {
public:
  int threads;
  uint64_t source;
  pf_execution_mode pf_mode;
  bool verify;
  int delta;
  int k_val;
  uint64_t seed64;
  boost::rand48 synch_gen;
  boost::uniform_int<uint64_t>* prand_vertex;  
  agm_work_stats work_stats;
  runtime_stats rt_stats;

  
  agm_instance_params(int _threads,
                      uint64_t s,
                      pf_execution_mode _mode,
                      bool _v,
                      int _delta,
                      int _k,
                      uint64_t _seed) : threads(_threads),
                                        source(s),
                                        pf_mode(_mode),
                                        verify(_v),
                                        delta(_delta),
                                        k_val(_k),
                                        seed64(_seed),
                                        prand_vertex(NULL),
                                        work_stats(threads),
                                        rt_stats(threads){
    synch_gen.seed(seed64);
  }

  std::string get_pf_execution_mode() {
    if (pf_mode == agm_pf_preorder)
      return "pre-order";
    else if (pf_mode == agm_pf_postorder)
      return "post-order";
    else if (pf_mode == agm_pf_splitted)
      return "splitted";
    else {
      error("Invalid processing function mode. Available : preorder, postorder, splitted");
      assert(false);
    }            
  }
  
  void print() {
    std::cout << "delta : " << delta << std::endl;
    std::cout << "k value : " << k_val << std::endl;    
    std::cout << "source : " << source << std::endl;
    std::cout << "verify : " << verify << std::endl;
    std::cout << "processing function execution mode: " << get_pf_execution_mode()
              << std::endl;
  }
  
  uint64_t get_source(uint64_t n) {
    if (source == undefined_source) {
      if (prand_vertex == NULL) {
	prand_vertex = new boost::uniform_int<uint64_t>(0, n-1);
      }
      // generate a random source
      return (*prand_vertex)(synch_gen);
    } else {
      return source;
    }
  }


  void print_summary() {
    work_stats.print_summary();    
  }

  std::string get_algorithm() {
    return "AGM-SSSP";
  }
};


class agm_params {

private:
  std::vector<agm_instance_params> params;
  int threads;
  std::vector<uint64_t> sources;
  pf_execution_mode pf_mode = agm_pf_splitted;  
  std::vector<int> k_vals;
  std::vector<int> deltas;    
  bool verify;
  uint64_t seed = 12345;
  
public:
  agm_params() : threads(1),
                 verify(false){}

  bool parse(int argc, char* argv[]){
        
    for (int i=0; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--threads") {
	threads = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--seed") {
	seed = boost::lexical_cast<uint64_t>( argv[i+1]);
      }
      
      if (arg == "--source") {
	sources = extract_params<uint64_t>( argv[i+1] );
      }      

      if (arg == "--pf-exec-mode") {
        std::string execmode = argv[i+1];
        if (execmode == "preorder")
          pf_mode = agm_pf_preorder;
        else if (execmode == "postorder")
          pf_mode = agm_pf_postorder;
        else if (execmode == "splitted")
          pf_mode = agm_pf_splitted;
        else {
          error("Invalid processing function execution mode. Available (preorder, postorder, splitted)");
          assert(false);
        }
      }            
      
      if (arg == "--k-val") {
	k_vals = extract_params<int>( argv[i+1] );
      }      

      if (arg == "--delta") {
	deltas = extract_params<int>( argv[i+1] );
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

      if (k_vals.empty())
        k_vals.push_back(2);

      if(deltas.empty())
        deltas.push_back(10);
      
      if (threads == -1) {
	std::cout << "[ERROR] Number of threads must be specified using --threads ..." << std::endl;
	assert(false);
      }

      BOOST_FOREACH (uint64_t s, sources) {
        BOOST_FOREACH (int k, k_vals) {
          BOOST_FOREACH (int d, deltas) {        
            agm_instance_params ainst(threads, s, pf_mode, verify, d, k, seed);
            params.push_back(ainst);
          }
        }
      }

    }
    
    return params;
  }
};

class AGMSSSPBaseExecutor {

public:
  template <typename Graph,
            typename EAGMConfig,
            typename MessageGenerator,
	    typename graph_create_params>
  time_type execute_sssp_eagm(Graph& g,
                              EAGMConfig& config,
                              amplusplus::transport& trans, 
                              MessageGenerator& msg_gen,
                              graph_create_params& gparams,
                              instance_params& runtime_params,
                              agm_instance_params& agm_params) {


    info("Running with following EAGM configuration;");
    config.print();

    boost::graph::agm::sssp_family<Graph> ssspes;
    
    // weights
    typedef typename boost::property_map<Graph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = boost::get(&WeightedEdge::weight, g);

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
    time_type exec_time = -1;
    if (agm_params.source != undefined_source) {
      exec_time = ssspes.execute_eagm(g,
                                      weight,
                                      rtgen,
                                      config,
                                      agm_params,
                                      agm_params.source,
                                      runtime_params,
                                      agm_params.work_stats,
                                      agm_params.verify);
    
      agm_params.work_stats.reduce_stats(exec_time);
      trans.set_nthreads(1);
      clear_thread_core_data();
      return exec_time;
      
    } else {
    
      while(exec_time == -1) {
        auto current_source = boost::vertex(agm_params.get_source(gparams.n), g);
        info("Executing source : ", current_source);
      
        exec_time = ssspes.execute_eagm(g,
                                        weight,
                                        rtgen,
                                        config,
                                        agm_params,
                                        current_source,
                                        runtime_params,
                                        agm_params.work_stats,
                                        agm_params.verify);
    
        agm_params.work_stats.reduce_stats(exec_time);
        trans.set_nthreads(1);
        clear_thread_core_data();
      }

      return exec_time;
    }    
  }


  template <typename Graph,
            typename EAGMConfig,
            typename MessageGenerator,
	    typename graph_create_params>
  time_type execute_sssp_kla_eagm(Graph& g,
                                  EAGMConfig& config,
                                  amplusplus::transport& trans, 
                                  MessageGenerator& msg_gen,
                                  graph_create_params& gparams,
                                  instance_params& runtime_params,
                                  agm_instance_params& agm_params) {


    info("Running with following EAGM configuration;");
    config.print();

    boost::graph::agm::sssp_family<Graph> ssspes;
    
    // weights
    typedef typename boost::property_map<Graph, weight_type WeightedEdge::*>::type WeightMap;
    WeightMap weight = boost::get(&WeightedEdge::weight, g);

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


    time_type exec_time = -1;
    if (agm_params.source != undefined_source) {
      exec_time = ssspes.execute_klevel_eagm(g,
                                             weight,
                                             rtgen,
                                             config,
                                             agm_params,
                                             agm_params.source,
                                             runtime_params,
                                             agm_params.work_stats,
                                             agm_params.verify);

    
      agm_params.work_stats.reduce_stats(exec_time);
      trans.set_nthreads(1);
      clear_thread_core_data();
      return exec_time;
      
    } else {
    
      while(exec_time == -1) {
        auto current_source = boost::vertex(agm_params.get_source(gparams.n), g);
        info("Executing source : ", current_source);
        exec_time = ssspes.execute_klevel_eagm(g,
                                               weight,
                                               rtgen,
                                               config,
                                               agm_params,
                                               current_source,
                                               runtime_params,
                                               agm_params.work_stats,
                                               agm_params.verify);

        
        agm_params.work_stats.reduce_stats(exec_time);
        trans.set_nthreads(1);
        clear_thread_core_data();
      }

      return exec_time;
    }
  }  
};


