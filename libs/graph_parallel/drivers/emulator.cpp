// Copyright (C) 2017 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//================== Gizmo AGM Emulator =====================//
// A performance emulator based on Abstract Graph Machine model.
//===========================================================//
#include "detail/bfs_executor.hpp"
#include "detail/emulator_executor.hpp"


template<typename graph_algorithm, 
	 typename algorithm_params, 
	 typename algo_instance_param>
class emulator {

public:
  emulator() {}

  int emulate(int argc, char* argv[]) {
    emulator_params eparams;
    eparams.parse(argc, argv);
    eparams.print();
    runtime_config_params rcp;
    rcp.parse(argc, argv);

    amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true,
							      1/*recvdepth*/,
							      1/*polls*/,
							      rcp.get_flow_control());
    amplusplus::transport trans = env.create_transport();

    _RANK = trans.rank();

    rcp.print();

    algorithm_params aparams;
    aparams.parse(argc, argv);
    
    assert(_RANK != -1);
    if (_RANK == 0) {
      std::cout << "============= Printing Algorithm Specific Parameters ================" << std::endl;
      aparams.print();
      std::cout << "======================================================================" << std::endl;
    }

    // Are we reading a graph from a file ?
    graph_reader_params readeparams;
    if (!readeparams.parse(argc, argv))
      return -1;

    // No we are not -- then generate the graph
    if (!readeparams.read_graph) {
      synthetic_generator sg(trans);

      graph_gen_params gparams;
      if (!gparams.parse(argc, argv))
	return -1;

      gparams.print();
      sg.generate_graph(gparams);
      
      if(_RANK == 0)
	std::cout << "[INFO] Done generating the graph ..." << std::endl;

      em_param_executor<graph_algorithm, 
		     algorithm_params, 
		     algo_instance_param, 
		     graph_gen_params> pexecutor;

      pexecutor.parameter_invoke(rcp, aparams, gparams, sg, env, trans, eparams);

    } else {
      graph_reader gr(trans, readeparams);
      // we are reading a graph from a file
      readeparams.print();
      gr.read_graph();

      em_param_executor<graph_algorithm, 
		     algorithm_params, 
		     algo_instance_param, 
		     graph_reader_params> pexecutor;

      pexecutor.parameter_invoke(rcp, aparams, readeparams, gr, env, trans, eparams);
    }

    return 0;

  }
};

int main(int argc, char* argv[]) {
  emulator<BFSExecutor, bfs_params, bfs_instance_params> bfsem;
  return bfsem.emulate(argc, argv);  
}
  
