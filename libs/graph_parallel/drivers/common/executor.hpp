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

#ifndef DRIVER_PBGL_EXECUTOR
#define DRIVER_PBGL_EXECUTOR
#include <iostream>
#include "synthetic_generator.hpp"
#include "graph_loader.hpp"
#include "parser.hpp"
#include <cstdlib>

template<typename graph_algorithm, 
	 typename algorithm_params, 
	 typename algo_instance_param,
	 typename graph_create_params>
class param_executor {

private:
  void print_input_summary(amplusplus::transport& trans,
			   graph_create_params& gparams,
			   instance_params& runtimeparams,
			   algo_instance_param& algoparam) {

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

    std::string coalescing_enabled = "Yes";
#ifdef NO_COALESCING
    coalescing_enabled = "No";
#endif    

    std::cout << "NODES : " << trans.size()
	      << ", THREADS : " << runtimeparams.threads
	      << ", GRAPH_TYPE : " << gparams.get_graph_type()
	      << gparams.get_generator_specific_input()
	      << ", ITERATIONS : " << runtimeparams.iterations
	      << ", COALESCING ENABLED : " << coalescing_enabled 
	      << ", COALESCING : " << runtimeparams.coalescing
	      << ", FLOW CONTROL : " << runtimeparams.flow_control
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
	      Graph& g,
	      MessageGenerator msg_gen,
	      graph_create_params& gparams,
	      instance_params& param,
	      algorithm_params& aparams) {
    // invoke the algorithm
    graph_algorithm algo;
    std::vector<algo_instance_param> alginstparams
      = aparams.get_instance_params();

    BOOST_FOREACH(algo_instance_param ainstp, alginstparams) {
      if (_RANK == 0) {
	std::cout << "==========================================" << std::endl;
	ainstp.print();
	std::cout << "==========================================" << std::endl;
      }
      // warming up
      for (int i=0; i < param.warm_up_iterations; ++i) {
	algo(g, trans, msg_gen, gparams, param, ainstp);
      }

      std::vector<time_type> all_times;
      for (int i=0; i < param.iterations; ++i) {
	time_type elapsed = algo(g, trans, msg_gen, gparams, param, ainstp);
	std::cout << "[INFO] Algorithm done in " << elapsed << "(sec) time." 
		  << std::endl;
	all_times.push_back(elapsed);
      }

      time_type mean = 0;
      time_type min = 0;
      time_type q1 = 0;
      time_type median = 0;
      time_type q3 = 0;
      time_type max = 0;
      time_type stddev = 0;

      time_statistics(all_times, 
		      mean, 
		      min,
		      q1,
		      median,
		      q3,
		      max,
		      stddev);

      if (_RANK == 0) {
	print_input_summary(trans, gparams, param, ainstp);
	std::cout << "MEAN : " << print_time(mean)
		  << " STDDEV : " << print_time(stddev)
		  << " MIN : " << print_time(min)
		  << " Q1 : " << print_time(q1)
		  << " MEDIAN : " << print_time(median)
		  << " Q3 : " << print_time(q3)
		  << " MAX : " << print_time(max)
		  << std::endl;
	ainstp.print_summary();
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
			amplusplus::transport& trans) {

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
		 gparams, param, aparams);

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
		 gparams, param, aparams);

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
		 gparams, param, aparams);

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
		 gparams, param, aparams);

	  break;
	}
      
	case rt_rook: {

	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::rook_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   amplusplus::rook_routing(trans.rank(), trans.size()));

	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams);

	  break;
	}

	case rt_hypercube: {
	  typedef amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> MessageGenerator;
	  MessageGenerator msg_gen(amplusplus::counter_coalesced_message_type_gen(param.coalescing), 
				   amplusplus::hypercube_routing(trans.rank(), trans.size()));

	  // invoke
	  invoke(trans, *(sg.graph()), msg_gen,
		 gparams, param, aparams);

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

template<typename graph_algorithm, 
	 typename algorithm_params, 
	 typename algo_instance_param>
class executor {

private:

public:
  int execute(int argc, char* argv[]) {
    std::cout << "Initializing AM++ environment ..." << std::endl;

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

      param_executor<graph_algorithm, 
		     algorithm_params, 
		     algo_instance_param, 
		     graph_gen_params> pexecutor;

      pexecutor.parameter_invoke(rcp, aparams, gparams, sg, env, trans);

    } else {
      graph_reader gr(trans, readeparams);
      // we are reading a graph from a file
      readeparams.print();
      gr.read_graph();

      param_executor<graph_algorithm, 
		     algorithm_params, 
		     algo_instance_param, 
		     graph_reader_params> pexecutor;

      pexecutor.parameter_invoke(rcp, aparams, readeparams, gr, env, trans);
    }
  }

};
#endif