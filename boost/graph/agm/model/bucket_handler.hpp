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

#ifndef __AGM_BUCKET_HANDLER__
#define __AGM_BUCKET_HANDLER__

#include <functional>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/model/eagm_buckets.hpp>
#include <boost/graph/agm/util/eagm_config.hpp>

#include <boost/parallel/append_buffer.hpp>

namespace boost { namespace graph { namespace agm {
      

template<typename EAGMConfig,
         typename Runtime,
         typename ProcessingFunction>
class bucket_handler {

  typedef typename Runtime::work_item_t work_item;

  typedef bucket_handler<EAGMConfig,
                         Runtime,
                         ProcessingFunction> self_type;


  typedef boost::graph::agm::eagm_traits<typename EAGMConfig::global_ordering_t,
                                         typename EAGMConfig::node_ordering_t,
                                         typename EAGMConfig::numa_ordering_t,
                                         typename EAGMConfig::thread_ordering_t,
                                         EAGMConfig,
                                         Runtime> eagm_traits_t;

  
  typedef boost::graph::agm::eagm_buckets<work_item,
                                          typename eagm_traits_t::root_level_bucket_trait>
  all_buckets_t;

  typedef append_buffer<work_item, 10u> InitialWorkItems;
  
public:
  bucket_handler(EAGMConfig& _config,
                 ProcessingFunction& _pf,
                 InitialWorkItems& _initial,
                 Runtime& _rt) :configs(_config),
                                pf(_pf),
                                initial(_initial),
                                runtime(_rt),
                                buckets(configs, runtime) {}

  void receive(const work_item& wi) {
    //    std::cout << "Invoking processing function ..................." << std::endl;
    const int tid = runtime.get_calling_thread_id();        
    pf(wi, tid, buckets);
  }
  
  void operator()(int tid) {
    typename EAGMConfig::global_ordering_t go;
    typename EAGMConfig::node_ordering_t no;
    typename EAGMConfig::numa_ordering_t nuo;
    typename EAGMConfig::thread_ordering_t to;

    // unfortunately we need to have this surrounded, TODO
    // think of a way to get rid of this
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      run(tid, go, no, nuo, to);
    }
  }

  time_type get_elapsed_time() {
    return elapsed_time;
  }

  template<typename T1,
           typename T2,
           typename T3,
           typename T4>
  void run(int tid, T1, T2, T3, T4) {

    runtime.initialize_per_thread(tid);

    runtime.synchronize();

    for (typename InitialWorkItems::size_type i = tid ;
         i < initial.size(); i+= runtime.get_nthreads()) {
      debug("pushing initial workitems");
      work_item& wi = initial[i];
      //pf(wi, tid, buckets);
      buckets.push(wi, tid);
    }
    
    runtime.synchronize();

    if (runtime.is_main_thread(tid)) {
      initial.clear();
    }
    runtime.wait_for_threads_to_reach_here(tid);

    // execution starts
    time_type start = get_time();    
    buckets.process_global_buckets(tid);
    time_type end = get_time();
    
    elapsed_time = (end-start);    
  }


  void run(int tid,
           CHAOTIC_ORDERING_T,
           CHAOTIC_ORDERING_T,
           CHAOTIC_ORDERING_T,
           CHAOTIC_ORDERING_T) {

    //    std::cout << "inside chaotic ordering" << std::endl;
    
    runtime.initialize_per_thread(tid);

    runtime.synchronize();

    // This is chaotic ordering, therefore, as soon
    // as we push an element execution starts    
    // execution starts
    time_type start = get_time();
    {
      runtime_epoch<Runtime> epoch(runtime);
      for (typename std::vector<work_item>::size_type i = tid ;
           i < initial.size(); i+= runtime.get_nthreads()) {
        work_item& wi = initial[i];
        buckets.push(wi, tid);
      }
    }
    time_type end = get_time();
    
    runtime.synchronize();
    
    elapsed_time = (end-start);    
  }

  void print_bucket_stats() {
    buckets.print_stats();
  }
  
  
private:
  EAGMConfig& configs;  
  ProcessingFunction& pf;
  InitialWorkItems& initial;  
  Runtime& runtime;
  all_buckets_t buckets;
  time_type elapsed_time;  
};


}}}
#endif
