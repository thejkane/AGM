#ifndef __AGM_PAGE_RANK_HPP__
#define __AGM_PAGE_RANK_HPP__

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <boost/graph/parallel/iteration_macros.hpp>
#include <boost/parallel/append_buffer.hpp>

#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/model/agm.hpp>
#include <boost/graph/agm/runtime/runtime.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>

namespace boost { namespace graph { namespace agm {
// definition of the CC work item set
typedef int Level;

template<typename Graph>
class kcore_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type KCoreType;

  // work item set definition
  // Every work item must have the time
  // destination, level
  typedef std::tuple<Vertex, KCoreType> WorkItem;

  
  // the processing function  
  template<typename State>
  struct kcore_pf {

  public:
    kcore_pf(const Graph& _rg,
                State& _st,
                agm_work_stats& _sr) : g(_rg),
                                       kcore(_st),
                                       stats(_sr){}

  private:
    const Graph& g;
    State& kcore;
    agm_work_stats& stats;

  public:
    template<typename buckets>
    void operator()(const WorkItem& wi,
                    int tid,
		    buckets& outset) {
      Vertex v = std::get<0>(wi);
      KCoreType core = std::get<1>(wi);
      
      KCoreType old_core = kcore[v], last_old_core;
      while(core < old_core) {
        last_old_core = old_core;
        old_core = boost::parallel::val_compare_and_swap(&kcore[v], old_core, core);
        if (last_old_core == old_core) {
#ifdef PBGL2_PRINT_WORK_STATS        
          if (old_core < boost::out_degree(v, g)) {
            stats.increment_invalidated(tid);
          } else
            stats.increment_useful(tid);            
#endif
          BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
            Vertex u = boost::target(e, g);
            WorkItem generated(u, core);
            outset.push(generated, tid);
          }

          return;          
        }
      }
#ifdef PBGL2_PRINT_WORK_STATS                
      stats.increment_rejected(tid);
#endif            
    }
  };


  template<typename append_buffer_t>
  struct threaded_initializer {
  private:
    int nthreads;
    Graph& graph;
    append_buffer_t* buffer;

  public:
    threaded_initializer(int _nthreads,
                         Graph& _g,
                         append_buffer_t* _buf) : nthreads(_nthreads),
                                                  graph(_g),
                                                  buffer(_buf) {}
    
    void operator()(int tid) {
      BGL_PARFORALL_VERTICES_T(u, graph, Graph, tid, nthreads) {
        BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
          Vertex v = target(e, graph);
          buffer->push_back(WorkItem(v, boost::out_degree(u, graph)));
        }
      }
    }
  };
  
public:
  template<typename KCoreMap>
  bool verify(Graph& g, KCoreMap& kcores) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    kcores.set_consistency_model(boost::parallel::cm_forward);
    kcores.set_max_ghost_cells(0);

    {
      amplusplus::scoped_epoch epoch(g.transport());
      
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(kcores, source(e, g));
	  get(kcores, target(e, g));
	}
      }
    }

    {	    
      amplusplus::scoped_epoch epoch(g.transport()); // at the moment get() sends a message
      BGL_FORALL_VERTICES_T(v, g, Graph) {
        KCoreType vk = get(kcores, v);
        // v must be connected to at least vk number of
        // vertices whos core is at least vk
        KCoreType i = 0;
	BGL_FORALL_ADJ_T(v, u, g, Graph) {
          if (get(kcores, u) >= vk) {
            ++i;
            if (i == vk)
              break;
          }          
	}

        if (i < vk) {
          std::cout << "[ERROR] Verification failed : Vertex : " << v
                    << ", has KCore = " << vk << ", but is not connected to neighbors who's"
                    << " kcore is at least " << vk << std::endl;
          std::cout << "Neighbors and their kcores : ";
          BGL_FORALL_ADJ_T(v, u, g, Graph) {
            std::cout << "{(" << v << ", " << u << ") : " << get(kcores, u) << "}, ";
          }

          assert(false);
        }
      }
    }
    
    kcores.clear(); // Clear memory used by ghost cells
    return true;
  }

  template<typename RuntimeModelGen,
           typename EAGMConfig,
           typename agm_instance_params>
  time_type execute_eagm(Graph& g,
                         RuntimeModelGen rtmodelgen,
                         EAGMConfig& config,
                         instance_params& runtime_params,
                         agm_instance_params& agm_params,
                         agm_work_stats& sr,
                         bool _verify=true) {
    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<KCoreType> kcoremap(num_vertices(g), 0);
    typedef boost::iterator_property_map<typename std::vector<KCoreType>::iterator, VertexIndexMap> KcoreMap;
    KcoreMap core_state(kcoremap.begin(), get(boost::vertex_index, g));

    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      put(core_state, v, boost::out_degree(v, g));
    }    
    
    info("Setting the initial work item set ...");
    // Initial work item set
    int nthreads = runtime_params.threads;
    typedef append_buffer<WorkItem, 10u> InitialBuffer;
    InitialBuffer initial;
    threaded_initializer<InitialBuffer> initializer(nthreads, g,
                                                    &initial); 
    boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
    for (int i = 0; i < nthreads - 1; ++i) {
      boost::thread thr(boost::ref(initializer), i + 1);
      threads[i].swap(thr);
    }
	  
    initializer(0);
    
    for (int i = 0; i < (nthreads - 1); ++i)
      threads[i].join();


    info("Setting the processing function ...");
    // Processing funcion
    typedef kcore_pf<KcoreMap> ProcessingFunction;
    ProcessingFunction pf(g, core_state, sr);
    
    // CC algorithm
    typedef eagm<Graph,
                 WorkItem,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen> kcore_eagm_t;

    kcore_eagm_t core(rtmodelgen,
                      config,
                      pf,
                      initial);

    
    info("Invoking Kcore algorithm ...");
    time_type elapsed = core(runtime_params);

    if (_verify) {
      verify(g, core_state);
    }

    return elapsed;
  }
};
}}}

#endif
