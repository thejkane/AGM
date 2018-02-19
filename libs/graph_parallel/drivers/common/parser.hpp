// Copyright (C) 2017 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine
#ifndef DRIVER_COMMANDLINE_PARSER
#define DRIVER_COMMANDLINE_PARSER
#include <vector>
#include <iostream>
#include <sstream>

#include <boost/graph/util/utils.hpp>


class instance_params {

private:

public:
  int threads;
  int iterations;
  int warm_up_iterations;
  bool enable_reductions;
  size_t reduction_cache;
  size_t coalescing;
  unsigned int poll_tasks;
  unsigned int flush;
  unsigned int eager;
  unsigned int recv_depth;
  routing_type routing;
  unsigned int flow_control;

  instance_params(int threads,
		  int iterations,
		  int warm,
		  bool redenabled,
		  size_t cache,
		  size_t coal,
		  unsigned int polls,
		  unsigned int flush,
		  unsigned int eager,
		  unsigned int recv,
		  routing_type rt,
		  unsigned int fc) : threads(threads),
				     iterations(iterations),
				     warm_up_iterations(warm),
				     enable_reductions(redenabled),
				     reduction_cache(cache),
				     coalescing(coal),
				     poll_tasks(polls),
				     flush(flush),
				     eager(eager),
    recv_depth(recv),
    routing(rt),
    flow_control(fc)
  {}

};

typedef std::vector<instance_params> common_parameters_t;

class runtime_config_params {

private:
  std::vector<int> thread_num_vals = {16};
  int iterations = 3;
  int warm_up_iterations = 1;
  bool enable_reductions = false;
  std::vector<size_t> reduction_cache_size = {10};
  std::vector<size_t> coalescing_size = {140000};
  std::vector<unsigned int> number_poll_task = {0};//disable
  std::vector<unsigned int> flush_freq = {20};
  std::vector<unsigned int> eager_limit = {3};
  std::vector<unsigned int> recv_depth = {0};//disable
  std::vector<routing_type> routing = {rt_none};
  unsigned int flow_control = 10;
  unsigned int reentry_count = 0;

  std::vector<instance_params> runtime_params;

public:
 
  std::vector<instance_params>& get_runtime_params() {
    return runtime_params;
  }

  unsigned int get_flow_control() {
    return flow_control;
  }

  unsigned int get_reentry_count() {
    return reentry_count;
  }


  template <typename vec>
  std::string print_list(vec param) {
    std::ostringstream oss;
    typename vec::iterator ite = param.begin();
    for (; ite != param.end(); ++ ite) {
      oss << (*ite) << ", ";
    }
    return oss.str();
  }

  void print() {
    assert(_RANK != -1);
    if (_RANK == 0) {
      std::cout << "==================== Printing Common Parameters =================" << std::endl;
      std::cout << "Threads = " << print_list(thread_num_vals) << std::endl;
      std::cout << "Iterations = " << iterations << std::endl;
      std::cout << "Warm up Iteration = " << warm_up_iterations << std::endl;
      std::cout << "Reductions enabled = " << enable_reductions << std::endl;
      std::cout << "Reduction Cache Sizes = " << print_list(reduction_cache_size) << std::endl;
      std::cout << "Coalescing Sizes = " << print_list(coalescing_size) << std::endl;
      std::cout << "Poll Tasks = " << print_list(number_poll_task) << std::endl;
      std::cout << "Flush Frequencies = " << print_list(flush_freq) << std::endl;
      std::cout << "Eager Limits = " << print_list(eager_limit) << std::endl;
      std::cout << "Receive Depths = " << print_list(recv_depth) << std::endl;
      std::cout << "Routing = " << print_list(routing) << std::endl;
      std::cout << "Flow Controls = " << flow_control << std::endl;
      std::cout << "Reentry Counts = " << reentry_count << std::endl;
      std::cout << "=================================================================" << std::endl;
    }
  }

  bool parse(int argc, char* argv[]) {
   
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if(arg == "--receive-depth") {
	recv_depth = extract_params<unsigned int>( argv[i+1] );
      }

      if (arg == "--threads") {
	thread_num_vals = extract_params<int>(argv[i+1]);
      }

      if (arg == "--iterations") {
	iterations = boost::lexical_cast<weight_type>( argv[i+1] );
      }

      if (arg == "--warm-up-iterations") {
	warm_up_iterations = boost::lexical_cast<weight_type>( argv[i+1] );
      }

      if (arg == "--coalescing-size") {
	coalescing_size = extract_params<size_t>( argv[i+1] );
      }

      if (arg == "--reduction-cache-size") {
	if (!enable_reductions) {
	  std::cerr << "Reductions must be enabled first. Use --enable-reductions to enable reductions"
		    << std::endl;
	  return false;
	} 

        reduction_cache_size = extract_params<size_t>( argv[i+1] );
      }

      if(arg == "--poll-task"){
	number_poll_task = extract_params<unsigned int>( argv[i+1] );
      }

      if(arg == "--flow-control"){
	flow_control = boost::lexical_cast<unsigned int>( argv[i+1] );
      }

      if(arg == "--reentry-count"){
	reentry_count = boost::lexical_cast<unsigned int>( argv[i+1] );
      }

      if (arg == "--flush") {
	flush_freq = extract_params<unsigned int> ( argv[i+1] );
      }

      if (arg == "--enable-reductions")
	enable_reductions = true;

      if (arg == "--routing") {
	std::string rt = boost::lexical_cast<std::string>( argv[i+1] );
	if (rt == "rt_none")
	  routing.push_back(rt_none);
	else if (rt == "rt_rook")
	  routing.push_back(rt_rook);
	else if (rt == "rt_hypercube")
	  routing.push_back(rt_hypercube);
	else 
	  std::cerr << "Invalid routing type. Available routing types -- rt_none, rt_rook and rt_hypercube"
		    << std::endl;
      }


    }

    if(coalescing_size.empty()) {
      std::cerr << "Provide a list of coalescing sizes: --coalescing-size x,y,..." << std::endl;
      return false;
    }

    if(recv_depth.empty()) {
      std::cerr << "Provide a list of receive depths: --receive-depth x,y,..." << std::endl;
      return false;
    }

    if(number_poll_task.empty()) {
      std::cerr << "Provide a list of poll tasks: --poll-task x,y,..." << std::endl;
      return false;
    }

    if(flush_freq.empty()) {
      std::cerr << "Provide a list of flush frequencies: --flush x,y,..." << std::endl;
      return false;
    }

    BOOST_FOREACH(routing_type rt, routing) {
      BOOST_FOREACH(int thread, thread_num_vals) {
	BOOST_FOREACH(size_t coalescing, coalescing_size) {
	  BOOST_FOREACH(unsigned int polls, number_poll_task) {
	    BOOST_FOREACH(unsigned int flush, flush_freq) {
	      BOOST_FOREACH(unsigned int eager, eager_limit) {
		BOOST_FOREACH(unsigned int recv, recv_depth) {
		  if (enable_reductions) {
		    BOOST_FOREACH(size_t reduction, reduction_cache_size) {
		      instance_params ip(thread,
					 iterations,
					 warm_up_iterations,
					 enable_reductions,
					 reduction,
					 coalescing,
					 polls,
					 flush,
					 eager,
					 recv,
					 rt,
					 flow_control);

		      runtime_params.push_back(ip);
		    }
		  } else {
		    instance_params ip(thread,
				       iterations,
				       warm_up_iterations,
				       enable_reductions,
				       0,
				       coalescing,
				       polls,
				       flush,
				       eager,
				       recv,
				       rt,
				       flow_control);

		    runtime_params.push_back(ip);
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    return true;
  }
};

class configuration_parser {

public:
  configuration_parser() {}

  void parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--scale") {
	genparams.scale = boost::lexical_cast<size_t>( argv[i+1] );
      }
    }
  }

private:
  graph_gen_params genparams;

};


#endif
