#ifndef __AGM_HPP__
#define __AGM_HPP__

#include <vector>
#include <tuple>
#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/agm/model/bucket_handler.hpp>

#define EMPTY_PF empty_processing_function

enum pf_execution_mode {
  agm_pf_preorder,
  agm_pf_postorder,
  agm_pf_splitted
};

namespace boost { namespace graph { namespace agm {

// An empty processing function.
// This will be used when post processing function or
// or the preprocessing function is not specified.
struct empty_processing_function {
public:
  template<typename work_item, typename out_set>
  void operator()(const work_item& wi, int tid, out_set& out) {
    out.push(wi, tid);
  }
};
      
template <typename Graph,
          typename WorkItem,
          typename PreOrderProcessingFunction,
          typename EAGMConfig,
	  typename RuntimeModelGen,
          typename PostOrderProcessingFunction>
class eagm {

  typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;

  class ReceiverHandler;
  // Machine model
  typedef typename RuntimeModelGen::template inner<WorkItem,
                                                   ReceiverHandler,
                                                   PostOrderProcessingFunction>::type RuntimeModel;
  // Bucket Handler
  typedef bucket_handler<EAGMConfig, RuntimeModel, PreOrderProcessingFunction> BucketHandler;
  // Initial Workitem Set
  typedef append_buffer<WorkItem, 10u> InitialWorkItems;

  class ReceiverHandler {
  private:
    BucketHandler* bucket_handler;
  public:
    ReceiverHandler() : bucket_handler(NULL){}
    ReceiverHandler(BucketHandler* _bh) : bucket_handler(_bh){}
    void operator() (const WorkItem& wi) const {
      bucket_handler->receive(wi);
    }
  };
  
public:
  // uses only pre-order processing function
  eagm (RuntimeModelGen& _rtgen,
        EAGMConfig& _config,
        PreOrderProcessingFunction& _pf,
        InitialWorkItems& initial_workitems) : empty_pf(),
                                               runtime(_rtgen, empty_pf),
                                               handler(_config, _pf,
                                                       initial_workitems, runtime),
                                               recv_handler(&handler){
    runtime.register_receiver(recv_handler);    
  }    
  // uses only post processing function
  eagm (RuntimeModelGen& _rtgen,
        EAGMConfig& _config,
        PostOrderProcessingFunction& _sendpf,
        InitialWorkItems& initial_workitems) : empty_pf(),
                                               runtime(_rtgen, _sendpf),
                                               handler(_config, empty_pf,
                                                       initial_workitems, runtime),
                                               recv_handler(&handler){
    runtime.register_receiver(recv_handler);
  }
  
  // uses both pre-order and post-order processing functions
  eagm (RuntimeModelGen& _rtgen,
        EAGMConfig& _config,
        PreOrderProcessingFunction& _pf,
        PostOrderProcessingFunction& _sendpf,
        InitialWorkItems& initial_workitems) : runtime(_rtgen, _sendpf),
                                               handler(_config, _pf,
                                                       initial_workitems, runtime),
                                               recv_handler(&handler){
    runtime.register_receiver(recv_handler);
  }

  time_type operator()(instance_params& runtime_params) {

    int nthreads = runtime_params.threads;

    runtime.set_transport_threads(nthreads);
    
    boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
    for (int i = 0; i < nthreads - 1; ++i) {
      boost::thread thr(boost::ref(handler), i + 1);
      threads[i].swap(thr);
    }
	  
    handler(0);
    
    for (int i = 0; i < (nthreads - 1); ++i)
      threads[i].join();

    runtime.set_transport_threads(1);

    return handler.get_elapsed_time();
  }

  void print_stats() {
    handler.print_bucket_stats();
  }

private:
  empty_processing_function empty_pf;
  RuntimeModel runtime;
  BucketHandler handler;
  ReceiverHandler recv_handler;  
};
      
      
template <typename Graph,
          typename WorkItem,
          typename PreOrderProcessingFunction,
          typename StrictWeakOrderingRelation,
	  typename RuntimeModelGen>
class agm {

  typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;

  // Machine model
  typedef typename RuntimeModelGen::template inner<WorkItem,
                                                   StrictWeakOrderingRelation,
                                                   PreOrderProcessingFunction>::type RuntimeModel;
public:
  agm (PreOrderProcessingFunction& _pf,
       StrictWeakOrderingRelation& _ord,
       RuntimeModelGen& _rtgen) : runtime(_pf, _ord, _rtgen) {}

  time_type operator()(std::vector<WorkItem>& initial_workitems,
           instance_params& runtime_params) {

    int nthreads = runtime_params.threads;

    runtime.set_threads(nthreads);
    
    boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
    for (int i = 0; i < nthreads - 1; ++i) {
      boost::thread thr(boost::ref(runtime), i + 1, initial_workitems);
      threads[i].swap(thr);
    }
	  
    runtime(0, initial_workitems);
    
    for (int i = 0; i < (nthreads - 1); ++i)
      threads[i].join();

    runtime.set_threads(1);

    return runtime.get_elapsed_time();
  }

private:
  RuntimeModel runtime;
};

}}}
#endif // __AGM_HPP
