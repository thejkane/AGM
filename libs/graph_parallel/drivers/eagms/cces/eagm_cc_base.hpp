// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

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

#include "../eagm_base.hpp"

#include <limits.h>
typedef boost::graph::agm::base_ordering base_ordering_t;

class agm_instance_params : public agm_instance_params_base {
public:
  int threads;
  id_distribution_t id_distribution;
  uint64_t delta;
  int kval;
  bool verify;
  agm_work_stats work_stats;
  runtime_stats rt_stats;

  
  agm_instance_params(pf_execution_mode _pfmode,
		      int _threads,
                      id_distribution_t& _idd,
                      uint64_t _delta,
		      int _k,
                      bool _v) : agm_instance_params_base(_pfmode),
				 threads(_threads),
                                 id_distribution(_idd),
                                 delta(_delta),
				 kval(_k),
                                 verify(_v),
                                 work_stats(threads),
                                 rt_stats(threads){}

  void print() {
    agm_instance_params_base::print();
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
    return "AGM-CC";
  }
};


class agm_params : public agm_params_base {

private:
  std::vector<agm_instance_params> params;
  int threads;
  std::vector<uint64_t> deltas;    
  id_distribution_t id_distribution = horizontal; // default  
  bool verify;
  std::vector<int> k_vals;
  
public:
  agm_params() : threads(1),
                 verify(false){}

  bool parse(int argc, char* argv[]){
    agm_params_base::parse(argc, argv);
        
    for (int i=0; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--threads") {
	threads = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--verify") {
        verify = true;
      }

      if (arg == "--k-val") {
	k_vals = extract_params<int>( argv[i+1] );
      }      

      if (arg == "--delta") {
	deltas = extract_params<uint64_t>( argv[i+1] );
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

      if (k_vals.empty())
        k_vals.push_back(1);

      if(deltas.empty())
        deltas.push_back(10);


      BOOST_FOREACH (int ks, k_vals) {
	BOOST_FOREACH (int delta, deltas) {        

	  agm_instance_params ainst(this->pf_mode,
				    threads, 
				    id_distribution, 
				    delta, 
				    ks, 
				    verify);
	  params.push_back(ainst);
	}
      }
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

    time_type exec_time = cces.execute_plain_eagm(g,
						  rtgen,
						  config,
						  runtime_params,
						  agm_params.work_stats,
						  iddist,
						  agm_params.verify);
    
    agm_params.work_stats.reduce_stats(exec_time);
    trans.set_nthreads(1);
    clear_thread_core_data();
    return exec_time;    
  }

  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename IdDistribution,
            typename EAGMConfig>
  time_type execute_with_level_eagm_config(Graph& g, 
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

    time_type exec_time = cces.execute_w_level_eagm(g,
						    rtgen,
						    config,
						    runtime_params,
						    agm_params.work_stats,
						    iddist,
						    agm_params.verify);
    
    agm_params.work_stats.reduce_stats(exec_time);
    trans.set_nthreads(1);
    clear_thread_core_data();
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


  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params,
            typename EAGMConfig>
  time_type select_id_distribution_w_level(Graph& g, 
					   amplusplus::transport& trans, 
					   MessageGenerator& msg_gen,
					   graph_create_params& gparams,
					   instance_params& runtime_params,
					   agm_instance_params& agm_params,
					   EAGMConfig& config) {

    if (agm_params.id_distribution == horizontal) {
      typedef row_id_distribution<Graph> id_dist_t;
      id_dist_t dist(g, trans.size());
      return execute_with_level_eagm_config(g,
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
      return execute_with_level_eagm_config(g,
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
