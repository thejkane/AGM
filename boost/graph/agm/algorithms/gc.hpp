#ifndef __AGM_GC_HPP
#define __AGM_GC_HPP

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

//===================== NOTE ===================================
// This is a hack! Need to fix it properly (TODO)
// This datatype should never be used but is needed to build the
// distributed property map
//===================== NOTE ===================================
namespace amplusplus {
  template <> struct make_mpi_datatype<std::vector< std::set<int> > > {MPI_Datatype get() const {return MPI_INT;}};

  template <> struct make_mpi_datatype<std::set<int>, void > {MPI_Datatype get() const {return MPI_INT;}};
}


namespace boost { namespace graph { namespace agm {
      
// definition of the GC work item set
typedef int Level;      

template<typename Graph,
         typename IdDistribution>
class gc_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;
  typedef int ColorType;

  // work item set definition
  // Every work item must have the time
  // destination, level
  typedef std::tuple<Vertex, Vertex, ColorType> WorkItem;
   
  // the processing function  
  template<typename State,
           typename NeighborCountMap,
           typename VertexNeighborColorMap>
  struct gc_pf {

  public:
    gc_pf(const Graph& _rg,
          State& _st,
          IdDistribution& _idd,
          NeighborCountMap& _lwerneighbr,
          NeighborCountMap& lwerfxedneighbrs,
          VertexNeighborColorMap& _vnbrclmp,
          agm_work_stats& _sr) : g(_rg),
                                 vcolor(_st),
                                 id_distribution(_idd),
                                 lower_neighbors(_lwerneighbr),
                                 lower_fixed_neighbors(lwerfxedneighbrs),
                                 vneighborcolors(_vnbrclmp),
                                 stats(_sr){}

  private:
    const Graph& g;
    State& vcolor;
    IdDistribution& id_distribution;
    NeighborCountMap& lower_neighbors;
    NeighborCountMap& lower_fixed_neighbors;
    VertexNeighborColorMap& vneighborcolors;
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
      ColorType source_color = std::get<2>(wi);

      assert(logical_id(source_vertex) < logical_id(dest_vertex));
      uint32_t expected = lower_fixed_neighbors[dest_vertex];
      uint32_t oldexpected = expected;
      uint32_t newval = oldexpected + 1;
        
      while(true) {
        newval = expected+1;
        if (__atomic_compare_exchange_n(&lower_fixed_neighbors[dest_vertex],
                                        &expected,
                                        newval,
                                        true,
                                        __ATOMIC_RELAXED,
                                        __ATOMIC_RELAXED)) {
          assert(newval == lower_fixed_neighbors[dest_vertex]);
          // update neighbor colors for this vertex
          // TODO this is not thread safe!
          vneighborcolors[dest_vertex].insert(source_color);
          if (newval == lower_neighbors[dest_vertex]) {
            vcolor[dest_vertex] = find_next_min_color(dest_vertex);
              
            std::set<Vertex> adjacencies;
            BGL_FORALL_OUTEDGES_T(dest_vertex, e, g, Graph) {
              Vertex u = target(e, g);
              if (logical_id(u) > logical_id(dest_vertex)) {
                if (adjacencies.insert(u).second) {
                  WorkItem wi(u, dest_vertex, vcolor[dest_vertex]);
                  outset.push(wi, tid);
                }
              }
            }
          }

          break;
        }                
      }
    }

    ColorType find_next_min_color(Vertex v) {
      ColorType min = 0;
      while(true) {
        auto ite = vneighborcolors[v].find(min);
        if (ite == vneighborcolors[v].end())
          return min;
        else
          ++min;
      }
    }
  };


public:
  template<typename GCStateMap>
  bool verify(Graph& g, GCStateMap& color, IdDistribution& idd) {
    amplusplus::transport& trans = g.transport();
    
    typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
    OwnerMap owner(get(boost::vertex_owner, g));

    color.set_consistency_model(boost::parallel::cm_forward || boost::parallel::cm_backward);
    color.set_max_ghost_cells(0);

    if (trans.rank() == 0) std::cout<<"Verifying GC results......";

    {
      amplusplus::scoped_epoch epoch(g.transport());
	  
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(color, target(e, g));
	}
      }
    }

    uint32_t scolor;
    std::set<Vertex> sources;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      // a neighbor of v cannot have the same color
      // as v
      bool source = true;
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
        Vertex u = target(e, g);
        if (idd(u) < idd(v))
          source = false;
        
        if (color[v] == color[u]) {
          std::cout << "[ERROR] Vertex v:" << v << ", and vertex u: " << u
                    << " have the same color. Vertex v color : " << color[v]
                    << ", and vertex u color : " << color[u]
                    << std::endl;
          assert(false);
        }
      }

      if (source) {
        scolor = color[v]; 
        sources.insert(v);
      }
    }

    // sources should have the same color value
    for (typename std::set<Vertex>::iterator ite = sources.begin();
         ite != sources.end(); ++ite) {
      assert(scolor == color[*ite]);
    }
  }


  template<typename append_buffer_t,
           typename NeighborCountMap,
           typename GCStateMap>
  struct threaded_initializer {
  private:
    int nthreads;
    Graph& graph;
    IdDistribution& idd;
    NeighborCountMap& lower_neighbors;
    GCStateMap& gc;
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
                         GCStateMap& _gc,
                         append_buffer_t* _buf) : nthreads(_nthreads),
                                                  graph(_g),
                                                  idd(_idd),
                                                  lower_neighbors(_lwerneighbrs),
                                                  gc(_gc),
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
          //TODO gc[u] = GC_FIX1;
          for (typename std::set<Vertex>::iterator ite = greater_adjacencies.begin();
               ite != greater_adjacencies.end(); ++ite) {
            buffer->push_back(WorkItem((*ite), u, gc[u]));
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
    std::vector<ColorType> gcmap(num_vertices(g), 0);
    typedef boost::iterator_property_map<typename std::vector<ColorType>::iterator, VertexIndexMap> GCStateMap;
    GCStateMap color_state(gcmap.begin(), get(boost::vertex_index, g));

    std::vector<uint32_t> lowerneighbrsmap(num_vertices(g), 0);
    std::vector<uint32_t> lowerfixedneighbrsmap(num_vertices(g), 0);    
    typedef boost::iterator_property_map<typename std::vector<uint32_t>::iterator, VertexIndexMap> NeighborCountMap;
    NeighborCountMap lower_neighbors(lowerneighbrsmap.begin(), get(boost::vertex_index, g));
    NeighborCountMap lower_fixed_neighbors(lowerfixedneighbrsmap.begin(), get(boost::vertex_index, g));    

    std::vector< std::set<ColorType> > vertexneighbormap(num_vertices(g));
    typedef boost::iterator_property_map<typename std::vector< std::set<ColorType> >::iterator, VertexIndexMap> VertexNeighborColorMap;
    VertexNeighborColorMap vneighborcolors(vertexneighbormap.begin(),
                                           get(boost::vertex_index, g));
    
    info("Setting the initial work item set ...");
    // Initial work item set
    int nthreads = runtime_params.threads;
    typedef append_buffer<WorkItem, 10u> InitialBuffer;
    InitialBuffer initial;
    threaded_initializer<InitialBuffer, NeighborCountMap, GCStateMap> initializer(nthreads, g,
                                                                                   iddist,
                                                                                   lower_neighbors,
                                                                                   color_state,
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
    typedef gc_pf<GCStateMap,
                  NeighborCountMap,
                  VertexNeighborColorMap> ProcessingFunction;
    ProcessingFunction pf(g,
                          color_state,
                          iddist,
                          lower_neighbors,
                          lower_fixed_neighbors,
                          vneighborcolors,
                          sr);
    
    // GC algorithm
    typedef eagm<Graph,
                 WorkItem,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen> gc_eagm_t;

    gc_eagm_t gcalgo(rtmodelgen,
                       config,
                       pf,
                       initial);

    
    info("Invoking GC algorithm ...");
    time_type elapsed = gcalgo(runtime_params);

    if (_verify) {
      verify(g, color_state, iddist);
    }

    return elapsed;
  }
  
};

}}}
#endif
