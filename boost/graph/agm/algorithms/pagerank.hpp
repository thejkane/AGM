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

#define DANGLING_VERTEX 0
#define NOT_DANGLING_VERTEX 1

#define ALL_DANGLING_RECV_TRUE 0
#define ALL_DANGLING_RECV_FALSE 1    
      
template<typename Graph,
         typename InDegreeMap>
class pagerank_family {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;
  typedef double PageRank;

  // work item set definition
  // Every work item must have the time
  // destination, level
  typedef int DangleVertexType;
  typedef std::tuple<Vertex, PageRank, DangleVertexType> WorkItem;
  typedef typename graph_traits<Graph>::degree_size_type degree_t;
  
  // the processing function  
  template<typename State,
           typename ReceiveEdgeCounterMap,
           typename DanglingRecvIndicatorMap,
           typename agm_instance_params>
  struct pagerank_pf {

  public:
    pagerank_pf(const Graph& _rg,
                InDegreeMap& _indeg,
                State& _st,
                State& _dynamicst,                
                State& _staticst,
                ReceiveEdgeCounterMap& _recvecountr,
                DanglingRecvIndicatorMap& _danglind,
                agm_instance_params& _insp,
                agm_work_stats& _sr) : g(_rg),
                                       indegrees(_indeg),
                                       pr_state(_st),
                                       pr_dynamic_state(_dynamicst),                                       
                                       pr_static_state(_staticst),
                                       recv_edge_counter(_recvecountr),
                                       dangling_recv_indicator(_danglind),
                                       params(_insp),
                                       stats(_sr){}

  private:
    const Graph& g;
    InDegreeMap& indegrees;
    State& pr_state;    
    State& pr_dynamic_state;
    State& pr_static_state;
    ReceiveEdgeCounterMap& recv_edge_counter;
    DanglingRecvIndicatorMap& dangling_recv_indicator;
    agm_instance_params& params;
    agm_work_stats& stats;

  public:
    template<typename buckets>
    void operator()(const WorkItem& wi,
                    int tid,
		    buckets& outset) {
      Vertex v = std::get<0>(wi);
      double inportion = std::get<1>(wi);
      DangleVertexType isdanglingrecv = std::get<2>(wi);

      // first atomically set the dangling indicator
      // this must be first action since if another thread
      // complete the receive count it may try to send the state as dangling
      if (isdanglingrecv == NOT_DANGLING_VERTEX) {
        // update dangling indicator to say that we have an incoming
        // message from a normal in vertex
        int expected = dangling_recv_indicator[v];
        int oldexpected = ALL_DANGLING_RECV_TRUE;
        int newstate = ALL_DANGLING_RECV_FALSE;

        while(oldexpected == expected) {
          if (__atomic_compare_exchange_n(&dangling_recv_indicator[v],
                                          &expected,
                                          newstate,
                                          true,
                                          __ATOMIC_RELAXED,
                                          __ATOMIC_RELAXED)) {
            assert(newstate == dangling_recv_indicator[v]);
          }
        }

        // increase dynamic (accumulating) page rank value
        // for dynamic portion, we will do the multiplication at the end
        // why this is not compiling in OSX ?
        auto newval = pr_dynamic_state[v]+inportion;
        auto expectedval = pr_dynamic_state[v];
        
        while(true) {
          expectedval = pr_dynamic_state[v];
          if (__atomic_compare_exchange(&pr_dynamic_state[v],
                                          &expectedval,
                                          &newval,
                                          true,
                                          __ATOMIC_RELAXED,
                                          __ATOMIC_RELAXED)) {
            break;
          }
        }

      } else {
        // rank update is from a dangling vertex,
        // update the static page rank value
        // this may lead to a race condition -- need to think!
        auto newval = pr_static_state[v]+inportion;
        auto expectedval = pr_static_state[v];
        
        while(true) {
          expectedval = pr_static_state[v];
          if (__atomic_compare_exchange(&pr_static_state[v],
                                        &expectedval,
                                        &newval,
                                        true,
                                        __ATOMIC_RELAXED,
                                        __ATOMIC_RELAXED)) {
            break;
          }
        }

      }
      
      // atomically increase the number of receives
      uint32_t receives_so_far = recv_edge_counter[v];
      
      while(true) {
        uint32_t newrecvs = receives_so_far+1;
        if (__atomic_compare_exchange_n(&recv_edge_counter[v],
                                        &receives_so_far,
                                        newrecvs,
                                        true,
                                        __ATOMIC_RELAXED,
                                        __ATOMIC_RELAXED)) {
          assert(newrecvs <= indegrees[v]);
          
          if (newrecvs == indegrees[v]) {

            // only one thread (the last thread) should enter this region
            // we have completed receives for the current iteration
            auto newrank = params.damping_factor * (pr_static_state[v] + pr_dynamic_state[v]);
            auto oldrank = pr_state[v];
            pr_state[v] = newrank;
            
            if ((newrank - oldrank) > params.epsillon) {
              // we are not done, we need to update out-edges
              // all the pageranks we receive are from danlging vertices ?
              // If yes, this vertex also, become a dangling vertex, then send
              // update neighbors saying current vertex is in dangling state
              auto portion = newrank / boost::out_degree(v, g);
              BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
                Vertex u = target(e, g);
                if (dangling_recv_indicator[v] == ALL_DANGLING_RECV_TRUE) {
                  // once the vertex goes to the dangling state, it will not receive
                  // any updates, therefore, the page rank value is final
                  WorkItem w(u, portion, DANGLING_VERTEX);
                  outset.push(w, tid);
                } else {
                  // if the vertex is not dangling, keep the static part
                  // reset the dynamic part
                  pr_dynamic_state[v] = 0;
                  // reset for next iteration
                  dangling_recv_indicator[v] = ALL_DANGLING_RECV_TRUE;
                  recv_edge_counter[v] = 0;
                  
                  WorkItem w(u, portion, NOT_DANGLING_VERTEX);
                  outset.push(w, tid);
                }
              }
            }
          }

          break;
        }
      }
    }
  };


  template<typename append_buffer_t>
  struct threaded_initializer {
  private:
    int nthreads;
    Graph& graph;
    InDegreeMap& indegrees;
    double initial_pr;
    append_buffer_t* buffer;

  public:
    threaded_initializer(int _nthreads,
                         Graph& _g,
                         InDegreeMap& _indeg,
                         double _init,
                         append_buffer_t* _buf) : nthreads(_nthreads),
                                                  graph(_g),
                                                  indegrees(_indeg),
                                                  initial_pr(_init),
                                                  buffer(_buf) {}
    
    void operator()(int tid) {
      BGL_PARFORALL_VERTICES_T(u, graph, Graph, tid, nthreads) {
        double portion = initial_pr / boost::out_degree(u, graph);
        BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
          Vertex v = target(e, graph);
          if (indegrees[u] == 0)
            buffer->push_back(WorkItem(v, portion, DANGLING_VERTEX));
          else
            buffer->push_back(WorkItem(v, portion, NOT_DANGLING_VERTEX));
        }
      }
    }
  };
  
public:
  template<typename PageRankState>
  bool verify(Graph& g, InDegreeMap& indegree, PageRankState& ranks) {
    // TODO
  }

  template<typename RuntimeModelGen,
           typename EAGMConfig,
           typename agm_instance_params>
  time_type execute_eagm(Graph& g,
                         InDegreeMap& indegrees,
                         RuntimeModelGen rtmodelgen,
                         EAGMConfig& config,
                         instance_params& runtime_params,
                         agm_instance_params& agm_params,
                         agm_work_stats& sr,
                         bool _verify=true) {
    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<PageRank> pagerankmap(num_vertices(g), 0);
    typedef boost::iterator_property_map<typename std::vector<PageRank>::iterator, VertexIndexMap> PageRankMap;
    PageRankMap pr_state(pagerankmap.begin(), get(boost::vertex_index, g));

    // PR dynamic state
    std::vector<PageRank> pagerankdynamicmap(num_vertices(g), 0);
    PageRankMap pr_dynamic_state(pagerankdynamicmap.begin(), get(boost::vertex_index, g));        

    // PR static state
    std::vector<PageRank> pagerankstaticmap(num_vertices(g), 0);
    PageRankMap pr_static_state(pagerankstaticmap.begin(), get(boost::vertex_index, g));    

    // To track receive in edges
    std::vector<uint32_t> recvedgecntermap(num_vertices(g), 0);
    typedef boost::iterator_property_map<typename std::vector<uint32_t>::iterator, VertexIndexMap> ReceiveEdgeCounterMap;
    ReceiveEdgeCounterMap recv_edge_counter(recvedgecntermap.begin(), get(boost::vertex_index, g));

    // To check whether all the incoming messages are from a dangling
    // vertices. If so the receiving vertex will go to dangling state from next
    // iteration
    std::vector<int> danglinrecv(num_vertices(g), ALL_DANGLING_RECV_TRUE);
    typedef boost::iterator_property_map<typename std::vector<int>::iterator, VertexIndexMap> DanglingRecvIndicatorMap;
    DanglingRecvIndicatorMap dangling_recv_ind(danglinrecv.begin(), get(boost::vertex_index, g));
    
    info("Setting the initial work item set ...");
    // Initial work item set
    int nthreads = runtime_params.threads;
    typedef append_buffer<WorkItem, 10u> InitialBuffer;
    InitialBuffer initial;
    threaded_initializer<InitialBuffer> initializer(nthreads,
                                                    g,
                                                    indegrees,
                                                    agm_params.initial_pr,
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
    typedef pagerank_pf<PageRankMap,
                        ReceiveEdgeCounterMap,
                        DanglingRecvIndicatorMap,
                        agm_instance_params> ProcessingFunction;
    
    ProcessingFunction pf(g,
                          indegrees,
                          pr_state,
                          pr_dynamic_state,
                          pr_static_state,
                          recv_edge_counter,
                          dangling_recv_ind,
                          agm_params,
                          sr);
    
    // CC algorithm
    typedef eagm<Graph,
                 WorkItem,
                 ProcessingFunction,
                 EAGMConfig,
                 RuntimeModelGen> pr_eagm_t;

    pr_eagm_t pralgo(rtmodelgen,
                     config,
                     pf,
                     initial);

    
    info("Invoking PageRank algorithm ...");
    time_type elapsed = pralgo(runtime_params);

    if (_verify) {
      verify(g, indegrees, pr_state);
    }

    return elapsed;
  }
};
}}}

#endif
