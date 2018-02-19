#include "common/synthetic_generator.hpp"
#include "common/graph_loader.hpp"
#include "common/parser.hpp"
#include "common/executor.hpp"
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
  int global_k;
  int node_k;
  int numa_k;
  int thread_k;
  uint64_t source;
  bool verify;
  agm_work_stats work_stats;
  base_ordering_t* global_ord;
  base_ordering_t* node_ord;
  base_ordering_t* numa_ord;
  base_ordering_t* thread_ord;
  runtime_stats rt_stats;

  
  agm_instance_params(int _threads,
                      int _gk,
                      int _nk,
                      int _nuk,
                      int _tk,                      
                      uint64_t s,
                      bool _v,
                      base_ordering_t* _globalo,
                      base_ordering_t* _nodeo,
                      base_ordering_t* _numao,
                      base_ordering_t* _threado) : threads(_threads),
                                                   global_k(_gk),
                                                   node_k(_nk),
                                                   numa_k(_nuk),
                                                   thread_k(_tk),
                                                   source(s),
                                                   verify(_v),
                                                   work_stats(threads),
                                                   global_ord(_globalo),
                                                   node_ord(_nodeo),
                                                   numa_ord(_numao),
                                                   thread_ord(_threado),
                                                   rt_stats(threads){}

  void print() {
    std::cout << "source : " << source << std::endl;
    std::cout << "global kla val : " << global_k << std::endl;
    std::cout << "node kla val : " << node_k << std::endl;
    std::cout << "numa kla val : " << numa_k << std::endl;
    std::cout << "thread kla val : " << thread_k << std::endl;    
    std::cout << "verify : " << verify << std::endl;    
  }

  void create_eagm_config() {
    
  }

  void print_summary() {
    work_stats.print_summary();    
  }

  std::string get_algorithm() {
    return "AGM";
  }
};


class agm_params {

private:
  std::vector<agm_instance_params> params;
  int threads;
  int global_k;
  int node_k;
  int numa_k;
  int thread_k;
  std::vector<uint64_t> sources;
  bool verify;
  base_ordering_t* global_ord;
  base_ordering_t* node_ord;
  base_ordering_t* numa_ord;
  base_ordering_t* thread_ord;
  
public:
  agm_params() : threads(1),
                 global_k(1),
                 node_k(1),
                 numa_k(1),
                 thread_k(1),
                 verify(false),
                 global_ord(NULL),
                 node_ord(NULL),
                 numa_ord(NULL),
                 thread_ord(NULL){}

  base_ordering_t* get_ordering(std::string arg) {

    base_ordering_t* pord
      = boost::graph::agm::ordering_type_map::get_ordering(arg);
    
    if (pord == NULL) {
      error("Invalid ordering. Available orderings :");
      boost::graph::agm::ordering_type_map::print();
      assert(false);
    } else {
      return pord;
    }
  }
  
  bool parse(int argc, char* argv[]){
        
    for (int i=0; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--threads") {
	threads = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--source") {
	sources = extract_params<uint64_t>( argv[i+1] );
      }      

      if (arg == "--global") {
        global_ord = get_ordering(argv[i+1]);
      }
      
      if (arg == "--node") {
        node_ord = get_ordering(argv[i+1]);
      }
      
      if (arg == "--numa") {
        numa_ord = get_ordering(argv[i+1]);
      }
      
      if (arg == "--thread") {
        thread_ord = get_ordering(argv[i+1]);
      }            
      
      if (arg == "--verify") {
        verify = true;
      }      

      if (arg == "--global-k") {
	global_k  = boost::lexical_cast<int>( argv[i+1] );
      }
      
      if (arg == "--node-k") {
	node_k  = boost::lexical_cast<int>( argv[i+1] );
      }
      
      if (arg == "--numa-k") {
	numa_k  = boost::lexical_cast<int>( argv[i+1] );
      }
      
      if (arg == "--thread-k") {
	thread_k  = boost::lexical_cast<int>( argv[i+1] );
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
  
      if (threads == -1) {
	std::cout << "[ERROR] Number of threads must be specified using --threads ..." << std::endl;
	assert(false);
      }

      BOOST_FOREACH (uint64_t s, sources) {

        agm_instance_params ainst(threads,
                                  global_k,
                                  node_k,
                                  numa_k,
                                  thread_k,
                                  s,
                                  verify,
                                  global_ord,
                                  node_ord,
                                  numa_ord,
                                  thread_ord);
        params.push_back(ainst);
      }

    }
    
    return params;
  }
};

class AGMBFSExecutor {

private:
  template <typename Graph,
            typename MessageGenerator,
            typename EAGMConfig,
	    typename graph_create_params>
  time_type invoke(const Graph& g, 
                   amplusplus::transport& trans, 
                   MessageGenerator& msg_gen,
                   EAGMConfig& config,
                   graph_create_params& gparams,
                   instance_params& runtime_params,
                   agm_instance_params& agm_params) {
    
    boost::graph::agm::bfs_family<Graph> bfses;

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

  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename GlobalOrdering,
            typename NodeOrdering,
            typename NumaOrdering>
  time_type select_thread_ordering (const Graph& g, 
                                    amplusplus::transport& trans, 
                                    MessageGenerator& msg_gen,
                                    graph_create_params& gparams,
                                    instance_params& runtime_params,
                                    agm_instance_params& agm_params,
                                    GlobalOrdering& globalord,
                                    NodeOrdering& nodeord,
                                    NumaOrdering& numaord) {
    
    base_ordering_t* thread = agm_params.thread_ord;
    assert(thread != NULL);
    
    if (thread->name() == CHAOTIC_ORDERING_T::ORDERING_NAME) {        
      CHAOTIC_ORDERING_T* threadord = static_cast<CHAOTIC_ORDERING_T*>(thread);
      auto config = boost::graph::agm::create_eagm_config(globalord,
                                                          nodeord,
                                                          numaord,
                                                          *threadord);
      invoke(g, trans,
             msg_gen,
             config,
             gparams,
             runtime_params,
             agm_params);
            
    } else if (thread->name() == LEVEL_ORDERING_T::ORDERING_NAME) {
      LEVEL_ORDERING_T* threadord = static_cast<LEVEL_ORDERING_T*>(thread);
      auto config = boost::graph::agm::create_eagm_config(globalord,
                                                          nodeord,
                                                          numaord,
                                                          *threadord);
      invoke(g, trans,
             msg_gen,
             config,
             gparams,
             runtime_params,
             agm_params);

    } else if (thread->name() == KLEVEL_ORDERING_T::ORDERING_NAME) {
      KLEVEL_ORDERING_T* threadord = static_cast<KLEVEL_ORDERING_T*>(thread);
      threadord->set_value(agm_params.thread_k);

      auto config = boost::graph::agm::create_eagm_config(globalord,
                                                          nodeord,
                                                          numaord,
                                                          *threadord);

      invoke(g, trans,
             msg_gen,
             config,
             gparams,
             runtime_params,
             agm_params);
            
    } else {
      error("Invalid thread ordering type : ", thread->name());
      assert(false);
    }
  }


  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename GlobalOrdering,
            typename NodeOrdering>
  time_type select_numa_ordering (const Graph& g, 
                                  amplusplus::transport& trans, 
                                  MessageGenerator& msg_gen,
                                  graph_create_params& gparams,
                                  instance_params& runtime_params,
                                  agm_instance_params& agm_params,
                                  GlobalOrdering& globalord,
                                  NodeOrdering& nodeord) {
    
    base_ordering_t* numa = agm_params.numa_ord;
    assert(numa != NULL);
    
    if (numa->name() == CHAOTIC_ORDERING_T::ORDERING_NAME) {        
      CHAOTIC_ORDERING_T* numaord = static_cast<CHAOTIC_ORDERING_T*>(numa);
      select_thread_ordering(g, trans,
                             msg_gen,
                             gparams,
                             runtime_params,
                             agm_params,
                             globalord,
                             nodeord,
                             *numaord);
               
    } else if (numa->name() == LEVEL_ORDERING_T::ORDERING_NAME) {
      LEVEL_ORDERING_T* numaord = static_cast<LEVEL_ORDERING_T*>(numa);
      select_thread_ordering(g, trans,
                             msg_gen,
                             gparams,
                             runtime_params,
                             agm_params,
                             globalord,
                             nodeord,
                             *numaord);
      
    } else if (numa->name() == KLEVEL_ORDERING_T::ORDERING_NAME) {
      KLEVEL_ORDERING_T* numaord = static_cast<KLEVEL_ORDERING_T*>(numa);
      select_thread_ordering(g, trans,
                             msg_gen,
                             gparams,
                             runtime_params,
                             agm_params,
                             globalord,
                             nodeord,
                             *numaord);
      
    } else {
      error("Invalid numa ordering type : ", numa->name());
      assert(false);
    }  
  }

  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename GlobalOrdering>
  time_type select_node_ordering (const Graph& g, 
                                  amplusplus::transport& trans, 
                                  MessageGenerator& msg_gen,
                                  graph_create_params& gparams,
                                  instance_params& runtime_params,
                                  agm_instance_params& agm_params,
                                  GlobalOrdering& globalord) {
    
    base_ordering_t* node = agm_params.node_ord;
    assert(node != NULL);
    
    if (node->name() == CHAOTIC_ORDERING_T::ORDERING_NAME) {        
      CHAOTIC_ORDERING_T* nodeord = static_cast<CHAOTIC_ORDERING_T*>(node);
      select_numa_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           globalord,
                           *nodeord);
        
    } else if (node->name() == LEVEL_ORDERING_T::ORDERING_NAME) {
      LEVEL_ORDERING_T* nodeord = static_cast<LEVEL_ORDERING_T*>(node);
      select_numa_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           globalord,
                           *nodeord);
      
    } else if (node->name() == KLEVEL_ORDERING_T::ORDERING_NAME) {
      KLEVEL_ORDERING_T* nodeord = static_cast<KLEVEL_ORDERING_T*>(node);
      select_numa_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           globalord,
                           *nodeord);
      
    } else {
      error("Invalid node ordering type : ", node->name());
      assert(false);
    }
  }  

  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params>
  time_type select_global_ordering (const Graph& g, 
                                  amplusplus::transport& trans, 
                                  MessageGenerator& msg_gen,
                                  graph_create_params& gparams,
                                  instance_params& runtime_params,
                                  agm_instance_params& agm_params) {

    base_ordering_t* global = agm_params.global_ord;
    assert(global != NULL);
    
    if (global->name() == CHAOTIC_ORDERING_T::ORDERING_NAME) {
      CHAOTIC_ORDERING_T* globalord = static_cast<CHAOTIC_ORDERING_T*>(global);
      select_node_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           *globalord);

    } else if (global->name() == DIJKSTRA_ORDERING_T::ORDERING_NAME) {
      DIJKSTRA_ORDERING_T* globalord = static_cast<DIJKSTRA_ORDERING_T*>(global);
      select_node_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           *globalord);
      
    } else if (global->name() == DIJKSTRA_ORDERING_STD_PQ_T::ORDERING_NAME) {
      DIJKSTRA_ORDERING_STD_PQ_T* globalord = static_cast<DIJKSTRA_ORDERING_STD_PQ_T*>(global);         select_node_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           *globalord);
         
    } else if (global->name() == DELTA_ORDERING_T::ORDERING_NAME) {
      DELTA_ORDERING_T* globalord = static_cast<DELTA_ORDERING_T*>(global);
      select_node_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           *globalord);
      
    } else if (global->name() == LEVEL_ORDERING_T::ORDERING_NAME) {
      LEVEL_ORDERING_T* globalord = static_cast<LEVEL_ORDERING_T*>(global);                             select_node_ordering(g, trans, msg_gen,
                           gparams,
                           runtime_params,
                           agm_params,
                           *globalord);
 
    } else {
      error("Invalid ordering type : ", global->name());
      assert(false);
    }

  }  
  
public:
  template <typename Graph, typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(const Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       agm_instance_params& agm_params) {

    CHAOTIC_ORDERING_T ch;
    select_numa_ordering(g, trans,
                         msg_gen,
                         gparams,
                         runtime_params,
                         agm_params,
                         ch,
                         ch);
  }
};

int main(int argc, char* argv[]) {

  boost::graph::agm::ordering_type_map::init();
  info("Starting AGM ...");
  std::cout << "[INFO] Printing core id for process ..." << std::endl;
  print_core_id();

  executor<AGMBFSExecutor, agm_params, agm_instance_params> agm_executor;
  agm_executor.execute(argc, argv);
  boost::graph::agm::ordering_type_map::clear();   
}
