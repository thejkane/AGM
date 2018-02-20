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

#ifndef __AGM_MIS_HPP
#define __AGM_MIS_HPP

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
      
#define MIS_UNFIX 0
#define MIS_FIX1 1
#define MIS_FIX0 2
      
// definition of the MIS work item set
typedef int Level;      

template<typename Graph,
         typename IdDistribution>
class mis_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;
  typedef int MISState;

  // work item set definition
  // Every work item must have the time
  // destination, level
  typedef std::tuple<Vertex, Vertex, MISState> WorkItem;
   
  // the processing function  
  template<typename State,
           typename NeighborCountMap>
  struct mis_pf {

  public:
    mis_pf(const Graph& _rg,
           State& _st,
           IdDistribution& _idd,
           NeighborCountMap& _lwerneighbr,
           NeighborCountMap& lwerfxedneighbrs,
           agm_work_stats& _sr) : g(_rg),
                                  vmis(_st),
                                  id_distribution(_idd),
                                  lower_neighbors(_lwerneighbr),
                                  lower_fixed_neighbors(lwerfxedneighbrs),
                                  stats(_sr){}

  private:
    const Graph& g;
    State& vmis;
    IdDistribution& id_distribution;
    NeighborCountMap& lower_neighbors;
    NeighborCountMap& lower_fixed_neighbors;
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
      
      Vertex dest_vertex = std::get<0>(wi);
      Vertex source_vertex = std::get<1>(wi);
      MISState source_state = std::get<2>(wi);
      
      if (vmis[dest_vertex] != MIS_UNFIX) {
        return;            
      }
      
      if (source_state == MIS_FIX1) {
        MISState expected = MIS_UNFIX;
        MISState newval = MIS_FIX0;

        while(expected == MIS_UNFIX) {
          if (__atomic_compare_exchange_n(&vmis[dest_vertex],
                                          &expected,
                                          newval,
                                          true,
                                          __ATOMIC_RELAXED,
                                          __ATOMIC_RELAXED)) {
            assert(vmis[dest_vertex] == MIS_FIX0);

            std::set<Vertex> adjacencies;
            BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
              Vertex u = target(e, g);
              if (logical_id(u) > logical_id(dest_vertex)) {
                if (adjacencies.insert(u).second) {
                  WorkItem wi(u, dest_vertex, vmis[dest_vertex]);
                  outset.push(wi, tid);
                }
              }
            }

            break;
          }
        }      
      } else if (source_state == MIS_FIX0) { // end of source_state == MIS_FIX1
        assert(logical_id(source_vertex) < logical_id(dest_vertex));
        uint32_t expected = lower_fixed_neighbors[dest_vertex];
        uint32_t oldexpected = expected;
        uint32_t newval = oldexpected + 1;
        
        while(true) {
          newval = expected + 1;
          if (__atomic_compare_exchange_n(&lower_fixed_neighbors[dest_vertex],
                                          &expected,
                                          newval,
                                          true,
                                          __ATOMIC_RELAXED,
                                          __ATOMIC_RELAXED)) {

            if (newval == lower_neighbors[dest_vertex]) {
              vmis[dest_vertex] = MIS_FIX1;
              
              std::set<Vertex> adjacencies;
              BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
                Vertex u = target(e, g);
                if (logical_id(u) > logical_id(dest_vertex)) {
                  if (adjacencies.insert(u).second) {
                    WorkItem wi(u, dest_vertex, vmis[dest_vertex]);
                    outset.push(wi, tid);
                  }
                }
              }
            }

            break;
          }                
        }        
      } else {
        assert(false);
      }
    }
  };


public:
  template<typename MISStateMap>
  bool verify(Graph& g, MISStateMap& mis) {
    amplusplus::transport& trans = g.transport();
    
    typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
    OwnerMap owner(get(boost::vertex_owner, g));

    mis.set_consistency_model(boost::parallel::cm_forward || boost::parallel::cm_backward);
    mis.set_max_ghost_cells(0);

    if (trans.rank() == 0) std::cout<<"Verifying MIS results......";

    {
      amplusplus::scoped_epoch epoch(g.transport());
	  
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(mis, target(e, g));
	}
      }
    }
	    
    bool result = true;
    int found_mis = 0;
#ifdef PRINT_DEBUG
    std::cout << "MIS = {";
#endif
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      if (get(mis, v) ==  MIS_FIX1) {// v in mis, none of the neigbours should be in mis
#ifdef PRINT_DEBUG
	if (v == 35) {
	  std::cout << "Printing all neighbors of " << v
		    << " -- [";

	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    std::cout << target(e, g) << ", ";
	  }
	  std::cout << "]" << std::endl;
	}

	std::cout << v << "-[";
#endif
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
#ifdef PRINT_DEBUG
	  std::cout << target(e, g) << ", ";
#endif
	  if (v == target(e, g))
	    continue;

	  if (get(mis, target(e, g)) == MIS_FIX1) {
	    std::cout << "[FAIL] Vertex : " << v << " and its neigbour " << target(e, g)
		      << " are both in MIS." 
		      << " Vertex value : " << get(mis, v)
		      << " neighbor value : " << get(mis, target(e, g)) 
		      << std::endl;
	    std::cout << "Neighbors of " << v << "- {";
	    BGL_FORALL_ADJ_T(v, u1, g, Graph) {
	      std::cout << "(" << u1 << ", " << get(mis, u1) << ")";
	    }
	    std::cout << "}" << std::endl;

	    auto k = target(e, g);
	    if (get(owner, k) == g.transport().rank()) {
	      std::cout << "Neighbors of " << k << "- {";
	      BGL_FORALL_ADJ_T(k, u2, g, Graph) {
		std::cout << "(" << u2 << ", " << get(mis, u2) << ")";
	      }
	      std::cout << "}" << std::endl;
	    }

	    result = false;
	  } else {
#ifdef PRINT_DEBUG
	    if (get(mis, target(e, g)) != MIS_FIX0) {
	      std::cout << "vertex : " 
			<< target(e, g) 
			<< ", is in " 
			<< get(mis, target(e, g))
			<< std::endl;
	    }
#endif

	    // cannot be MIS_UNFIX
	    assert(get(mis, target(e, g)) == MIS_FIX0);
	    found_mis = 1;
	  }
	}
#ifdef PRINT_DEBUG
	std::cout << "], ";
#endif
      } else {
	// if v is not in mis, at least one of its neighbors must
	// be in mis
	if (get(mis, v) != MIS_FIX0) {
	  std::cout << "[ERROR] Vertex - " << v << " is in " << get(mis, v)
		    << std::endl;
	}

	if (get(mis, v) != MIS_FIX0) {
	  std::cout << "Error : " << v << " is not in MIS_FIX0. But in " << get(mis, v)
		    << ", owner of " << v << ":" << get(owner, v)
		    << ", current rank : " << _RANK << std::endl;
	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    auto k = target(e, g);
	    std::cout << "neigbor : " << k << " mis : " << get(mis, k)
		      << std::endl;
	  }
	}
	assert(get(mis, v) == MIS_FIX0);

	bool inmis = false;
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  if (v == target(e, g))
	    continue;

	  if (get(mis, target(e, g)) == MIS_FIX1) {
	    inmis = true;
	    break;
	  }
	}
      
	if (!inmis) {
	  std::cout << "Vertex : " << v << " and none of its neighbors is in MIS" << std::endl;
	  assert(false);
	}
      }
    }
#ifdef PRINT_DEBUG
    std::cout << "}" << std::endl;
#endif

    // Run a reduction on found mis.
    // We cannot expect every rank will have found_mis true.
    // E.g:- rank0 - {0,1}, rank1 - {281474976710656, 281474976710657}
    // If every vertex is connect to each other, a vertex (Lets say 0) will
    // be in one rank. In the second rank all vertices are out of mis.
    int or_found_mis = 0;
    MPI_Allreduce(&found_mis, &or_found_mis, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (or_found_mis == 0) {
      std::cout << "Did not find any MIS" <<std::endl;
      result = false;
    }

    return result;
    
    
  }


  template<typename append_buffer_t,
           typename NeighborCountMap,
           typename MISStateMap>
  struct threaded_initializer {
  private:
    int nthreads;
    Graph& graph;
    IdDistribution& idd;
    NeighborCountMap& lower_neighbors;
    MISStateMap& mis;
    append_buffer_t* buffer;

    template<typename SizeType>
    SizeType logical_id(SizeType k) {
      return idd(k);
    }
        
  public:
    threaded_initializer(int _nthreads,
                         Graph& _g,
                         IdDistribution& _idd,
                         NeighborCountMap& _lwerneighbrs,
                         MISStateMap& _mis,
                         append_buffer_t* _buf) : nthreads(_nthreads),
                                                  graph(_g),
                                                  idd(_idd),
                                                  lower_neighbors(_lwerneighbrs),
                                                  mis(_mis),
                                                  buffer(_buf) {}
    
    void operator()(int tid) {
      BGL_PARFORALL_VERTICES_T(u, graph, Graph, tid, nthreads) {
        std::set<Vertex> adjacencies;
        std::set<Vertex> greater_adjacencies;        
        uint32_t count = 0;
        BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
          Vertex v = target(e, graph);
          if (logical_id(v) < logical_id(u)) {
            if (adjacencies.insert(v).second) { //to avoid parallel edges
              ++count;
            }
          } else {
            greater_adjacencies.insert(v);
          }
        }

        if (count == 0) {
          mis[u] = MIS_FIX1;
          for (typename std::set<Vertex>::iterator ite = greater_adjacencies.begin();
               ite != greater_adjacencies.end(); ++ite) {
            buffer->push_back(WorkItem((*ite), u, mis[u]));
          }
        }

        lower_neighbors[u] = count;
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
    std::vector<MISState> mismap(num_vertices(g), MIS_UNFIX);
    typedef boost::iterator_property_map<typename std::vector<MISState>::iterator, VertexIndexMap> MISStateMap;
    MISStateMap mis_state(mismap.begin(), get(boost::vertex_index, g));

    std::vector<uint32_t> lowerneighbrsmap(num_vertices(g), 0);
    std::vector<uint32_t> lowerfixedneighbrsmap(num_vertices(g), 0);    
    typedef boost::iterator_property_map<typename std::vector<uint32_t>::iterator, VertexIndexMap> NeighborCountMap;
    NeighborCountMap lower_neighbors(lowerneighbrsmap.begin(), get(boost::vertex_index, g));
    NeighborCountMap lower_fixed_neighbors(lowerfixedneighbrsmap.begin(), get(boost::vertex_index, g));    
    
    
    info("Setting the initial work item set ...");
    // Initial work item set
    int nthreads = runtime_params.threads;
    typedef append_buffer<WorkItem, 10u> InitialBuffer;
    InitialBuffer initial;
    threaded_initializer<InitialBuffer, NeighborCountMap, MISStateMap> initializer(nthreads, g,
                                                                                   iddist,
                                                                                   lower_neighbors,
                                                                                   mis_state,
                                                                                   &initial);
    
    boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
    for (int i = 0; i < nthreads - 1; ++i) {
      boost::thread thr(boost::ref(initializer), i + 1);
      threads[i].swap(thr);
    }
	  
    initializer(0);
    
    for (int i = 0; i < (nthreads - 1); ++i)
      threads[i].join();

    assert(initial.size() != 0);    

    info("Setting the processing function ...");
    // Processing funcion
    typedef mis_pf<MISStateMap, NeighborCountMap> ProcessingFunction;
    ProcessingFunction pf(g,
                          mis_state,
                          iddist,
                          lower_neighbors,
                          lower_fixed_neighbors,
                          sr);
    
    // MIS algorithm
    typedef eagm<Graph,
                 WorkItem,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen> mis_eagm_t;

    mis_eagm_t misalgo(rtmodelgen,
                       config,
                       pf,
                       initial);

    
    info("Invoking MIS algorithm ...");
    time_type elapsed = misalgo(runtime_params);

    if (_verify) {
      verify(g, mis_state);
    }

    return elapsed;
  }
  
};

}}}
#endif
