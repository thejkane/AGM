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

#ifndef __AGM_BFS_HPP
#define __AGM_BFS_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/model/agm.hpp>
#include <boost/graph/agm/runtime/runtime.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>

namespace boost { namespace graph { namespace agm {
      
// definition of the BFS work item set
typedef int Level;

template<typename Graph>
class bfs_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;

  // work item set definition
  // Every work item must have the time
  // destination, level
  typedef std::tuple<Vertex, Level> WorkItem;

  // the processing function  
  template<typename State>
  struct bfs_pf {

  public:
    bfs_pf(const Graph& _rg, State& _st, agm_work_stats& _sr) : g(_rg),
						       vlevel(_st),
						       stats(_sr){}

  private:
    const Graph& g;
    State& vlevel;
    agm_work_stats& stats;

  public:

    template<typename buckets>
    void operator()(const WorkItem& wi,
                    int tid,
		    buckets& outset) {
      
      Vertex v = std::get<0>(wi);

      std::cout << "============================================================" << std::endl;
      info("vertex : " , v);
      int level = std::get<1>(wi);
      int old_level = vlevel[v], last_old_level;

      info("new level : ", level);
      debug("old level : " , old_level);
      
      while(level < old_level) {
        last_old_level = old_level;
        old_level = boost::parallel::val_compare_and_swap(&vlevel[v], old_level, level);

        if (last_old_level == old_level) {
          
#ifdef PBGL2_PRINT_WORK_STATS        
          if (old_level < INT_MAX) {
            stats.increment_invalidated(tid);
          } else
            stats.increment_useful(tid);            

#endif
          BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
            Vertex u = boost::target(e, g);
            WorkItem generated(u, (level+1));
            outset.push(generated, tid);
            debug("pushing to outset : ", u);
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
  template<typename DistMap>
  bool verify(const Graph& g, Vertex s, DistMap& dists) {


    dists.set_consistency_model(boost::parallel::cm_forward);
    dists.set_max_ghost_cells(0);

    if (g.transport().rank() == 0)
      info("Verifying distances ....");      

    {
      amplusplus::scoped_epoch epoch(g.transport());
	  
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
          boost::get(dists, target(e, g));
	}
      }
    }

    
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	if (dists[source(e, g)] == INT_MAX) {
	  if (dists[target(e, g)] != INT_MAX) {
	    std::cout << "[ERROR] Source unvisited, therefore target must also be unvisited. But Target -- (" << target(e, g) << ", " << dists[target(e, g)] << ") -> Source -- ("
		      << source(e, g) << ", " << dists[source(e, g)] << ") " << std::endl;
	    std::abort();
	  }
	} else if (dists[target(e, g)] > (dists[source(e, g)] + 1)) {
	  std::cout << "[ERROR] Target distance is greater than source distance + 1: Target -- (" << target(e, g) << ", " << dists[target(e, g)] << ") -> Source -- ("
		    << source(e, g) << ", " << dists[source(e, g)] << ") " << std::endl;	  
	  std::abort();
	}
      }
    }

    info("Verification successful ...");
    return true;
  }


  template<typename RuntimeModelGen,
           typename EAGMConfig>
  time_type execute_eagm(const Graph& g,
                         RuntimeModelGen rtmodelgen,
                         EAGMConfig& config,
                         instance_params& runtime_params,
                         agm_work_stats& sr, 
                         bool _verify=true) {

    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<int> distmap(num_vertices(g), INT_MAX);
    typedef boost::iterator_property_map<typename std::vector<int>::iterator, VertexIndexMap> DistMap;
    DistMap distance_state(distmap.begin(), get(boost::vertex_index, g));

    // source
    Vertex source = 2;

    info("Setting the initial work item set ...");
    // Initial work item set
    std::vector<WorkItem> initial;
    initial.push_back(WorkItem(source, 0));

    info("Setting the processing function ...");
    // Processing funcion
    typedef bfs_pf<DistMap> ProcessingFunction;
    ProcessingFunction pf(g, distance_state, sr);
    
    // BFS algorithm
    typedef eagm<Graph,
                 WorkItem,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen> bfs_eagm_t;

    bfs_eagm_t bfsalgo(rtmodelgen,
                       config,
                       pf,
                       initial);

    
    info("Invoking BFS algorithm ...");
    time_type elapsed = bfsalgo(runtime_params);

    if (_verify) {
      verify(g, source, distance_state);
    }

    return elapsed;
  }
  
  template<typename RuntimeModelGen,
           typename StrictWeakOrdering>
  time_type execute(const Graph& g,
                    RuntimeModelGen rtmodelgen,
                    StrictWeakOrdering ordering,
		    instance_params& runtime_params,
		    agm_work_stats& sr, 
		    bool _verify=true) {

    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<int> distmap(num_vertices(g), INT_MAX);
    typedef boost::iterator_property_map<typename std::vector<int>::iterator, VertexIndexMap> DistMap;
    DistMap distance_state(distmap.begin(), get(boost::vertex_index, g));

    // source
    Vertex source = 2;

    info("Setting the initial work item set ...");
    // Initial work item set
    std::vector<WorkItem> initial;
    initial.push_back(WorkItem(source, 0));

    info("Setting the processing function ...");
    // Processing funcion
    typedef bfs_pf<DistMap> ProcessingFunction;
    ProcessingFunction pf(g, distance_state, sr);
    
    // BFS algorithm
    typedef agm<Graph,
                WorkItem,
                ProcessingFunction,
                StrictWeakOrdering,
                RuntimeModelGen> bfs_agm_t;

    bfs_agm_t bfsalgo(pf,
                      ordering,
                      rtmodelgen);
    
    info("Invoking BFS algorithm ...");
    time_type elapsed = bfsalgo(initial, runtime_params);

    if (_verify) {
      verify(g, source, distance_state);
    }

    return elapsed;
  }
};

}}}
#endif
