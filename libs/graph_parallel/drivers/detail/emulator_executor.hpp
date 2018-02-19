#ifndef EMULATOR_PBGL_EXECUTOR
#define EMULATOR_PBGL_EXECUTOR
#include <iostream>
#include "../common/synthetic_generator.hpp"
#include "../common/graph_loader.hpp"
#include "../common/parser.hpp"
#include "stat_reader.hpp"
#include <cstdlib>

enum em_graph_type {
  erdos,
  rmat1,
  rmat2,
  road,
  social
};


class emulator_params {

public:
  int abort_level;
  em_graph_type gtype;
  std::string get_graph_type() {
    if (gtype == erdos)
      return "Erdos-Renyi Graph";
    else if (gtype == rmat1)
      return "RMAT-1";
    else if (gtype == rmat2)
      return "RMAT-2";
    else if (gtype == road)
      return "Road";
    else if (gtype == social)
      return "Social";
    else
      return "Undefined";
  }

public:
  emulator_params():abort_level(std::numeric_limits<int>::max()),
		    gtype(erdos){}

  void parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--em-erdos")
	gtype = erdos;

      if (arg == "--em-abort-level")
	abort_level = boost::lexical_cast<int>( argv[i+1]);
    }
  }

  void print() {
    std::cout << "=============== Printing emulator params ===========================" << std::endl;
    std::cout << "Abort level : " << abort_level << std::endl;
    std::cout << "Graph Type : " << get_graph_type() << std::endl;
    std::cout << "====================================================================" << std::endl;
  }
};

template<typename graph_algorithm, 
	 typename algorithm_params, 
	 typename algo_instance_param,
	 typename graph_create_params>
class em_param_executor {

private:
  void print_input_summary(amplusplus::transport& trans,
			   graph_create_params& gparams,
			   instance_params& runtimeparams,
			   algo_instance_param& algoparam,
			   emulator_params& emparams) {

#ifdef DISABLE_SELF_SEND_CHECK
    bool self_send_disable = true;
#else
    bool self_send_disable = false;
#endif

    const char* threadlevel = std::getenv("MPICH_MAX_THREAD_SAFETY");
    if (threadlevel == NULL)
      threadlevel = "not defined";

    const char* asyncprogress = std::getenv("MPICH_NEMESIS_ASYNC_PROGRESS");
    if (asyncprogress == NULL)
      asyncprogress = "not defined";



    std::cout << "NODES : " << trans.size()
	      << ", THREADS : " << runtimeparams.threads
	      << ", GRAPH_TYPE : " << gparams.get_graph_type()
	      << gparams.get_generator_specific_input()
	      << ", ITERATIONS : " << runtimeparams.iterations
	      << ", COALESCING : " << runtimeparams.coalescing
	      << ", FLUSH : " << runtimeparams.flush
	      << ", REDUCTION ENABLED : " << runtimeparams.enable_reductions
	      << ", REDUCTION_CACHE : " << runtimeparams.reduction_cache
	      << ", SELF_SEND_DISABLE : " << self_send_disable
	      << ", MPI_THREAD_LEVEL : " << threadlevel
	      << ", MPI_ASYNC_PROGRESS : " << asyncprogress
	      << ", ALGORITHM : " << algoparam.get_algorithm()
	      << std::endl;

  }

  template<typename Graph, typename MessageGenerator>
  void invoke(amplusplus::transport& trans,
	      const Graph& g,
	      MessageGenerator msg_gen,
	      graph_create_params& gparams,
	      instance_params& param,
	      algorithm_params& aparams,
	      emulator_params& emparams) {
    // invoke the algorithm
    graph_algorithm algo;
    std::vector<algo_instance_param> alginstparams
      = aparams.get_instance_params();

    BOOST_FOREACH(algo_instance_param ainstp, alginstparams) {
      std::cout << "==========================================" << std::endl;
      ainstp.print();
      std::cout << "==========================================" << std::endl;

      typedef stat_reader<graph_create_params,
			  instance_params,
			  algorithm_params,
			  emulator_params> StatReader_t;
      StatReader_t sreader(gparams, param, aparams, emparams);

      // emulation
      // For erdos graphs we need at least two rising readings and two decending readins
      WorkStatData_t usefuldiff;
      WorkStatData_t invalidatediff;
      WorkStatData_t rejectediff;


      int64_t cuseful = 0;
      int64_t cinvalid = 0;
      int64_t crejected = 0;

      int level = 0;

      while(true) {
	algo.set_abort_level(level);
	time_type elapsed = algo(g, trans, msg_gen, gparams, param, ainstp);

	std::cout << "cuseful=" << cuseful << "algo.get_useful_work=" << algo.get_useful_work() << std::endl;
	usefuldiff.push_back(algo.get_useful_work() - cuseful);
	cuseful = algo.get_useful_work();

	invalidatediff.push_back(algo.get_invalidated_work() - cinvalid);
	cinvalid = algo.get_invalidated_work();

	rejectediff.push_back(algo.get_rejected_work() - crejected);
	crejected = algo.get_rejected_work();

	algo.reset_work_stats();
	std::cout << "[INFO] Emulation step :" << level << " done in " << elapsed << "(sec) time." 
		  << std::endl;

	if (emparams.abort_level == level)
	  break;
	else
	  ++level;
      }

      // process statistics
      sreader.process(usefuldiff,
		      invalidatediff,
		      rejectediff);

      if (_RANK == 0) {
	print_input_summary(trans, gparams, param, ainstp, emparams);
      }

    }
  }

public:
  template<typename graph_creator>
  void parameter_invoke(runtime_config_params& rcp,
			algorithm_params& aparams,
			graph_create_params& gparams,
			graph_creator& sg,
			amplusplus::environment& env,
			amplusplus::transport& trans,
			emulator_params& emparams) {

    common_parameters_t& rparams = rcp.get_runtime_params();
    common_parameters_t::iterator ite = rparams.begin();
    for(; ite != rparams.end(); ++ite) {
      instance_params& param = (*ite);

      if (param.poll_tasks != 0) {
	env.template downcast_to_impl<amplusplus::detail::mpi_environment_obj>()->set_poll_tasks((*ite).poll_tasks);
      }

      if (param.recv_depth != 0) {
	env.template downcast_to_impl<amplusplus::detail::mpi_environment_obj>()->set_recv_depth((*ite).recv_depth);
      }

      if (param.enable_reductions) {
	switch(param.routing) {
	case rt_none: {
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, 
							 amplusplus::no_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   param.reduction_cache, 
				   amplusplus::no_routing(trans.rank(), trans.size()));
	 
	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams, emparams);

	  break;
	}
      
	case rt_rook: {
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, 
							 amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   param.reduction_cache, 
				   amplusplus::rook_routing(trans.rank(), trans.size()));

	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams, emparams);

	  break;
	}

	case rt_hypercube: {
	  typedef amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, 
							 amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   param.reduction_cache, 
				   amplusplus::hypercube_routing(trans.rank(), trans.size()));
	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams, emparams);

	  break;
	}
	
	}
      } else {
	// reductions are disabled
	switch(param.routing) {
	case rt_none: {
	  typedef amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> MessageGenerator;
	  MessageGenerator msg_gen((amplusplus::counter_coalesced_message_type_gen(param.coalescing)));
	 
	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams, emparams);

	  break;
	}
      
	case rt_rook: {

	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   amplusplus::rook_routing(trans.rank(), trans.size()));

	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams, emparams);

	  break;
	}

	case rt_hypercube: {
	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   amplusplus::hypercube_routing(trans.rank(), trans.size()));

	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams, emparams);

	  break;
	}
	
	}
      }
    }

    assert(rparams.size() >= 1);
    sg.finalize();
    amplusplus::clear_mpi_datatype_map();
  }
};

#endif
