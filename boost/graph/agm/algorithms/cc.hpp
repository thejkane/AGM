#ifndef __AGM_CC_HPP
#define __AGM_CC_HPP

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

template<typename Graph,
         typename IdDistribution>
class cc_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;
  typedef Vertex Component;

  // work item set definition
  // Every work item must have the time
  // destination, level
  typedef std::tuple<Vertex, Component> WorkItem;
   
  // the processing function  
  template<typename State>
  struct cc_pf {

  public:
    cc_pf(const Graph& _rg,
          State& _st,
          IdDistribution& _idd,
          agm_work_stats& _sr) : g(_rg),
                                 vcomponent(_st),
                                 id_distribution(_idd),
                                 stats(_sr){}

  private:
    const Graph& g;
    State& vcomponent;
    IdDistribution& id_distribution;
    agm_work_stats& stats;

    template<typename SizeType>
    SizeType logical_id(SizeType k) {
      return id_distribution(k);
    }
    
    
  public:
    template<typename buckets>
    void operator()(const WorkItem& wi,
                    int tid,
		    buckets& outset) {
      
      Vertex v = std::get<0>(wi);
      int component = std::get<1>(wi);
      int old_component = vcomponent[v], last_old_component;
      
      while(component < old_component) {
        last_old_component = old_component;
        old_component = boost::parallel::val_compare_and_swap(&vcomponent[v], old_component, component);

        if (last_old_component == old_component) {
          
#ifdef PBGL2_PRINT_WORK_STATS        
          if (old_component < logical_id(v)) {
            stats.increment_invalidated(tid);
          } else
            stats.increment_useful(tid);            

#endif
          std::set<Vertex> adjacencies;
          bool haslowernbr = false;
          BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
            Vertex u = boost::target(e, g);
            if (logical_id(u) > component) { // ignore self-loops
              adjacencies.insert(u);
            } else if(logical_id(u) < component) {
              // v has a neighbor that is lower than v_component. Therefore, stop the search from v_component
              haslowernbr = true;
              break;
            }
          }

          if (!haslowernbr) {
            typename std::set<Vertex>::iterator ite = adjacencies.begin();
            for(; ite != adjacencies.end(); ++ite) {
              WorkItem generated((*ite), component);
              outset.push(generated, tid);              
            }
          }

          return;
        }
      }
      
#ifdef PBGL2_PRINT_WORK_STATS                
      stats.increment_rejected(tid);
#endif      
    }
  };


public:
  template<typename ComponentMap>
  bool verify(Graph& g, ComponentMap& components) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    components.set_consistency_model(boost::parallel::cm_forward);
    components.set_max_ghost_cells(0);
	      
    {
      amplusplus::scoped_epoch epoch(g.transport());
      
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(components, source(e, g));
	  get(components, target(e, g));
	}
      }
    }

    {	    
      amplusplus::scoped_epoch epoch(g.transport()); // at the moment get() sends a message
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_ADJ_T(v, u, g, Graph) {
	  //	  std::cout << "verifying vertex v : " << v << std::endl;
	  //#ifdef PRINT_DEBUG
	  if (get(components, v) != get(components, u)) 
	    std::cout << "Component of " << v << " : " << get(components, v)
		      << " component of " << u << " : " << get(components, u)
		      << std::endl;

	  assert(get(components, v) == get(components, u)); 
	}
      }
    }
    
    components.clear(); // Clear memory used by ghost cells
    return true;
  }


  template<typename append_buffer_t>
  struct threaded_initializer {
  private:
    int nthreads;
    Graph& graph;
    IdDistribution& idd;
    append_buffer_t* buffer;

    template<typename SizeType>
    SizeType logical_id(SizeType k) {
      return idd(k);
    }
        
  public:
    threaded_initializer(int _nthreads,
                         Graph& _g,
                         IdDistribution& _idd,
                         append_buffer_t* _buf) : nthreads(_nthreads),
                                                  graph(_g),
                                                  idd(_idd),
                                                  buffer(_buf) {}
    
    void operator()(int tid) {

      BGL_PARFORALL_VERTICES_T(u, graph, Graph, tid, nthreads) {
        Vertex min_neighbor = logical_id(u);
        std::set<Vertex> adjacencies;
        BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
          Vertex v = target(e, graph);
          if (u != v) { // ignore self loops
            if (adjacencies.insert(v).second) {
              // check whether u is the minimum out of its neighbors
              if (logical_id(v) < min_neighbor) {
                min_neighbor = logical_id(v);
                break;
              }
            }
          }
        }

        if (min_neighbor == logical_id(u)) {
          auto comp = logical_id(u);
          for (typename std::set<Vertex>::iterator ite = adjacencies.begin();
               ite != adjacencies.end(); ++ite) {
            assert(comp < logical_id(*ite));
            buffer->push_back(WorkItem((*ite), comp));
          }
        }
      }
    }
  };
  
  
  template<typename RuntimeModelGen,
           typename EAGMConfig>
  time_type execute_eagm(Graph& g,
                         RuntimeModelGen rtmodelgen,
                         EAGMConfig& config,
                         instance_params& runtime_params,
                         agm_work_stats& sr,
                         IdDistribution& iddist,
                         bool _verify=true) {

    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<Component> componentmap(num_vertices(g), std::numeric_limits<Component>::max());
    typedef boost::iterator_property_map<typename std::vector<Component>::iterator, VertexIndexMap> ComponentMap;
    ComponentMap component_state(componentmap.begin(), get(boost::vertex_index, g));
    
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      put(component_state, v, iddist(v));
    }

    info("Setting the initial work item set ...");
    // Initial work item set
    int nthreads = runtime_params.threads;
    typedef append_buffer<WorkItem, 10u> InitialBuffer;
    InitialBuffer initial;
    threaded_initializer<InitialBuffer> initializer(nthreads, g, iddist, &initial); 
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
    typedef cc_pf<ComponentMap> ProcessingFunction;
    ProcessingFunction pf(g, component_state, iddist, sr);
    
    // CC algorithm
    typedef eagm<Graph,
                 WorkItem,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen> cc_eagm_t;

    cc_eagm_t ccalgo(rtmodelgen,
                     config,
                     pf,
                     initial);

    
    info("Invoking CC algorithm ...");
    time_type elapsed = ccalgo(runtime_params);

    if (_verify) {
      verify(g, component_state);
    }

    return elapsed;
  }
  
};

}}}
#endif
