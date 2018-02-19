// Copyright (C) 2017 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine
#ifndef DRIVER_COMMON_UTIL
#define DRIVER_COMMON_UTIL

#include <mutex>
#include <boost/lexical_cast.hpp>
#include <boost/graph/distributed/time_calc.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
#include <mpi.h>
#include <utmpx.h>
#include <unistd.h>
#ifdef __APPLE__
#include <sys/unistd.h>
#else
#include <asm/unistd.h>
#endif
#include <sys/syscall.h>

#include <set>
#include <mpi.h>
#include <iomanip> 

int _RANK = -1;


#ifdef __APPLE__
int sched_getcpu() {
  return 0;
}
#endif

// testing functions
void print_core_id() {
#ifndef __APPLE__
  printf( "cpu = %d\n", sched_getcpu() );
#endif
}

static inline long my_gettid(void) { 
#ifndef __APPLE__
  return syscall(SYS_gettid); 
#else
  return 0;
#endif
}

void print_tid() {
#ifndef __APPLE__
  printf( "tid = %d\n", my_gettid() );
#endif
}

void print_tid_coreid() {
#ifndef __APPLE__
  printf( "tid = %d, coreid= %d\n", my_gettid(), sched_getcpu() );
#endif
}

int get_num_cores() {
#ifndef __APPLE__
  return sysconf(_SC_NPROCESSORS_ONLN);
#else
  return 1;
#endif
}

std::mutex set_mtx;
std::set<int> thread_core_data;
void validate_thread_core_relation() {
#ifndef __APPLE__
  int coreid = sched_getcpu();
  set_mtx.lock();
  if (!thread_core_data.insert(coreid).second) {
    std::cout << "core : " << coreid << " is used by two threads" << std::endl;
    set_mtx.unlock();
    assert(false);
  } else
    set_mtx.unlock();
#endif  
}


void clear_thread_core_data() {
  thread_core_data.clear();
}

int pin(int core) {
#ifndef __APPLE__
  int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
  if (core < 0 || core >= num_cores) {
    std::cerr << "Invalid core " << core << std::endl;
    return 1;
  }

  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core, &cpuset);
  // In Linux underneath boost is using pthread. May not compatible
  // with other platforms
  pthread_t current_thread = pthread_self();    
  return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
#else
  return 0;
#endif
}



template<typename T>
void time_statistics(std::vector<T>& data, 
		     T& m, 
		     T& minimum,
		     T& q1,
		     T& median,
		     T& q3,
		     T& maximum,
		     T& stddev) {

  if (data.empty())
    return;

  using namespace boost::accumulators;
  accumulator_set<T, stats<tag::variance> > acc;
  std::for_each(data.begin(), data.end(), boost::bind<void>(boost::ref(acc), _1));
  
  m = mean(acc);
  stddev = sqrt(variance(acc));

  auto const min = 0;
  auto const Q1 = data.size() / 4;
  auto const Q2 = data.size() / 2;
  auto const Q3 = Q1 + Q2;
  auto const max = data.size() - 1;


  std::nth_element(data.begin(), data.begin() + min, data.end());
  minimum = data[min];
  std::nth_element(data.begin(), data.begin() + Q1, data.end());
  q1 = data[Q1];
  std::nth_element(data.begin() + Q1 + 1, data.begin() + Q2, data.end());
  median = data[Q2];
  std::nth_element(data.begin() + Q2 + 1, data.begin() + Q3, data.end()); 
  q3 = data[Q3];
  std::nth_element(data.begin() + Q3 + 1, data.begin() + max, data.end()); 
  maximum = data[max];
}


template<typename Result>
std::vector<Result> extract_params(const std::string args) {
  size_t d = 0, d2;

  std::vector<Result> r;
  while ((d2 = args.find(',', d)) != std::string::npos) {
    r.push_back(boost::lexical_cast<Result>(args.substr(d, d2 - d)));
    d = d2 + 1;
  }

  r.push_back(boost::lexical_cast<Result>(args.substr(d, args.length())));

  return r;
}


// Needed so atomics are available on MacOS
// Edge weights
typedef int32_t weight_type; 
struct WeightedEdge {
  WeightedEdge(weight_type weight = 0) : weight(weight) { }

  weight_type weight;
};

// Distributions
enum distribution_t {
  block,
  cyclic
};

// Routings
enum routing_type {rt_none, rt_hypercube, rt_rook};


std::string print_time(time_type t)
{
  std::ostringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
  return out.str();
}

//=============== Simple Logger ================================//
void log(std::ostream& stream) {
  std::cout << std::endl;
}

template<typename H, typename... Ts>
void log(std::ostream& stream, H& h, Ts&... ts) {
  stream << h;
  log(stream, ts...);
}

template<typename... Ts>
void log(Ts... ts) {
  log(std::cout, ts...);
}

template<typename... Ts>
void info(Ts... ts) {
  if (_RANK == 0)
    log("[INFO] ", ts...);
}

template<typename... Ts>
void error(Ts... ts) {
  if (_RANK == 0)
    log("[ERROR] ", ts...);
}

template<typename... Ts>
void debug(Ts... ts) {
  if (_RANK == 0)
    log("[DEBUG] ", ts...);
}

template<typename... Ts>
void warn(Ts... ts) {
  if (_RANK == 0)
    log("[WARNING] ", ts...);
}

template<typename... Ts>
void all_rank_info(Ts... ts) {
  log("[INFO] ", ts...);
}

template<typename... Ts>
void all_rank_error(Ts... ts) {
  log("[ERROR] ", ts...);
}

template<typename... Ts>
void all_rank_debug(Ts... ts) {
  log("[DEBUG] ", ts...);
}

template<typename... Ts>
void all_rank_warn(Ts... ts) {
  log("[WARNING] ", ts...);
}

#endif
