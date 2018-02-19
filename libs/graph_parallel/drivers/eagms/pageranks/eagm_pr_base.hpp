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
#include <boost/graph/agm/algorithms/pagerank.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/runtime/ampp_runtime.hpp>

#include "in_degree_calculator.hpp"

#include <limits.h>
typedef boost::graph::agm::base_ordering base_ordering_t;

class agm_instance_params {
public:
  int threads;
  id_distribution_t id_distribution;
  bool verify;
  double damping_factor;
  double epsillon;
  int iterations;
  double initial_pr;
  agm_work_stats work_stats;
  runtime_stats rt_stats;

  
  agm_instance_params(int _threads,
                      id_distribution_t& _idd,
                      bool _v,
                      double _damp,
                      double _epslon,
                      int _iter,
                      double _init) : threads(_threads),
                                      id_distribution(_idd),
                                      verify(_v),
                                      damping_factor(_damp),
                                      epsillon(_epslon),
                                      iterations(_iter),
                                      initial_pr(_init),
                                      work_stats(threads),
                                      rt_stats(threads){}

  void print() {
    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;

    std::cout << "verify : " << verify << std::endl;    
  }


  void print_summary() {
    work_stats.print_summary();    
  }

  std::string get_algorithm() {
    return "AGM-PageRank-Chaotic";
  }
};


class agm_params {

private:
  std::vector<agm_instance_params> params;
  int threads;
  double damping_factor = 0.15;
  double epsillon = 0.2;
  int iterations = 2;
  double initial_pr = 1;
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

      if (arg == "--damping-factor") {
        damping_factor = boost::lexical_cast<double>( argv[i+1]);
      }

      if (arg == "--epsillon") {
        epsillon = boost::lexical_cast<double>( argv[i+1]);
      }
      
      if (arg == "--initial-pagerank") {
        initial_pr = boost::lexical_cast<double>( argv[i+1]);
      }
      
      if (arg == "--iterations") {
        iterations = boost::lexical_cast<int>( argv[i+1]);
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

      agm_instance_params ainst(threads,
                                id_distribution,
                                verify,
                                damping_factor,
                                epsillon,
                                iterations,
                                initial_pr);
      params.push_back(ainst);
    }
    
    return params;
  }
};



class AGMPRBaseExecutor {
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

    // we using a csr graph, therefore, we do not have
    // access to vertex in degrees. We calculate vertex
    // in degrees as a preprpocessing step.
    // TODO This should really be moved to graph generation
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    typedef std::vector< std::atomic<uint32_t> > AtomicInDegreeMap;
    AtomicInDegreeMap inedgeatimics(num_vertices(g));

    std::vector< uint32_t > inedgecounts(num_vertices(g));
    typedef boost::iterator_property_map<typename std::vector< uint32_t >::iterator, VertexIndexMap> InEdgeCounterMap;
    InEdgeCounterMap in_edge_counter_map(inedgecounts.begin(), get(boost::vertex_index, g));

    typedef VertexInDegreeCalculator<InEdgeCounterMap,
                                     AtomicInDegreeMap,
                                     Graph> VertexDegreeCalculatorType;

    int nthreads = runtime_params.threads;
    trans.set_nthreads(nthreads);
    VertexDegreeCalculatorType calc(g,
                                    trans,
                                    in_edge_counter_map,
                                    inedgeatimics,
                                    nthreads);
    
    boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
    for (int i = 0; i < nthreads - 1; ++i) {
      boost::thread thr(boost::ref(calc), i + 1);
      threads[i].swap(thr);
    }
	  
    calc(0);
    
    for (int i = 0; i < (nthreads - 1); ++i)
      threads[i].join();

    inedgeatimics.clear();
    trans.set_nthreads(1);

    /*    BGL_FORALL_VERTICES_T(v, g, Graph) {
      assert(boost::out_degree(v, g) == in_edge_counter_map[v]);
      }*/
    
    boost::graph::agm::pagerank_family<Graph, InEdgeCounterMap> prs;    
    
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

    time_type exec_time = prs.execute_eagm(g,
                                           in_edge_counter_map,
                                           rtgen,
                                           config,
                                           runtime_params,
                                           agm_params,
                                           agm_params.work_stats,
                                           agm_params.verify);
    
    agm_params.work_stats.reduce_stats(exec_time);
    return exec_time;    
  }

  virtual ~AGMPRBaseExecutor(){}
  
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
