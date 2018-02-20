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