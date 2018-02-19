#include "eagm_sssp_base.hpp"

class AGMSSSPDeltaExecutor : public AGMSSSPBaseExecutor {

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
    DELTA_ORDERING_T delta(agm_params.delta);
    auto config = boost::graph::agm::create_eagm_config(delta,
                                                        ch,
                                                        ch,
                                                        ch);

    return execute_sssp_eagm(g,
                             config,
                             trans,
                             msg_gen,
                             gparams,
                             runtime_params,
                             agm_params);
    
  }
};

int main(int argc, char* argv[]) {

  boost::graph::agm::ordering_type_map::init();
  info("Starting AGM ...");
  std::cout << "[INFO] Printing core id for process ..." << std::endl;
  print_core_id();

  executor<AGMSSSPDeltaExecutor, agm_params, agm_instance_params> agm_executor;
  agm_executor.execute(argc, argv);
  boost::graph::agm::ordering_type_map::clear();   
}
