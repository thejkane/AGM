#ifndef __IN_DEGREE_CALC__
#define __IN_DEGREE_CALC__
template<typename InDegreeMap,
         typename AtomicInDegreeMap,
         typename Graph,
         typename MessageGenerator =
         amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
struct VertexInDegreeCalculator {
  
  typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;

  typedef typename boost::graph_traits<Graph>::degree_size_type degree_sz_t;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  //typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
  
  //typedef boost::iterator_property_map<typename std::vector< std::atomic<uint32_t> >::iterator, VertexIndexMap> AtomicInDegreeMap;
  
  typedef Vertex in_edge_data_t; 

  struct vertex_degree_handler {
  public:
    vertex_degree_handler() : pg(NULL),
                              in_degrees(NULL) {}
    vertex_degree_handler(Graph* _pg,
                          AtomicInDegreeMap* _deg) : pg(_pg),
                                                     in_degrees(_deg){}
    void operator() (const in_edge_data_t& v) const {
      (*in_degrees)[pg->distribution().local(v)].fetch_add(1);
    }    
  private:
    Graph* pg;
    AtomicInDegreeMap* in_degrees;
  };
  
  
  typedef typename MessageGenerator::template
  call_result<in_edge_data_t,
              vertex_degree_handler, 
              OwnerMap,
              amplusplus::no_reduction_t >::type
  RelaxMessage;

  typedef VertexInDegreeCalculator<InDegreeMap, Graph, MessageGenerator> self_type;
  
public:
  VertexInDegreeCalculator(Graph& _g,
                           amplusplus::transport& _trans,
                           InDegreeMap& _map,
                           AtomicInDegreeMap& _atomicmap,
                           int _threads,
                           MessageGenerator msg_gen=MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12))) :
    dummy((amplusplus::register_mpi_datatype<in_edge_data_t>(), 0)),
    g(_g),
    transport(_trans),
    owner(get(boost::vertex_owner, g)),
    in_degree(_map),
    atomic_degrees(_atomicmap),
    nthreads(_threads),
    relax_msg(msg_gen, transport, owner, amplusplus::no_reduction) {

    relax_msg.set_handler(vertex_degree_handler(&g,
                                                &atomic_degrees));
  }

public:
  void operator()(int tid) {
    AMPLUSPLUS_WITH_THREAD_ID(tid) {    
      {
        amplusplus::scoped_epoch epoch(transport);
        BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
          BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
            Vertex v = target(e, g);
            relax_msg.send(v);
          }
        }
      }

      // load counter
      BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
        in_degree[u] = atomic_degrees[g.distribution().local(u)].load();
      }
    }
  }
  
private:
  int dummy;
  Graph g;
  amplusplus::transport& transport;  
  const OwnerMap& owner;
  InDegreeMap& in_degree;
  AtomicInDegreeMap& atomic_degrees;
  int nthreads;
  RelaxMessage relax_msg;
  
};

#endif
