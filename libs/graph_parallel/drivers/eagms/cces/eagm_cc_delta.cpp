#include "eagm_cc_base.hpp"

// Static initializations
using general_orderings = boost::graph::agm::ordering_type_map;
general_orderings::ordering_map_t general_orderings::ordering_map = general_orderings::create_ordering_map();

class AGMCCDeltaExecutor : public AGMCCBaseExecutor {

public:  
  template <typename Graph,
            typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       agm_instance_params& agm_params) {
    
    CHAOTIC_ORDERING_T ch;
    DELTA_ORDERING_T deltaord(agm_params.delta);
    auto config = boost::graph::agm::create_eagm_config(ch,
                                                        deltaord,
                                                        ch,
                                                        ch);    
    return select_id_distribution(g,
				  trans,
				  msg_gen,
				  gparams,
				  runtime_params,
				  agm_params,
                                  config);
  }
};

int main(int argc, char* argv[]) {

  boost::graph::agm::ordering_type_map::init();
  info("Starting AGM ...");
  std::cout << "[INFO] Printing core id for process ..." << std::endl;
  print_core_id();

  executor<AGMCCDeltaExecutor, agm_params, agm_instance_params> agm_executor;
  agm_executor.execute(argc, argv);
  boost::graph::agm::ordering_type_map::clear();   
}
