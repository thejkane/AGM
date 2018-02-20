// Copyright 2010-2013 The Trustees of Indiana University.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#include <config.h>

#include <boost/config.hpp>
#define IS_MPI_TRANSPORT_mpi 1
#define IS_MPI_TRANSPORT_gasnet 0
#define IS_MPI_TRANSPORT_shm 0
#define IS_MPI_TRANSPORT BOOST_JOIN(IS_MPI_TRANSPORT_, TRANSPORT)
#define IS_SHM_TRANSPORT_mpi 0
#define IS_SHM_TRANSPORT_gasnet 0
#define IS_SHM_TRANSPORT_shm 1
#define IS_SHM_TRANSPORT BOOST_JOIN(IS_SHM_TRANSPORT_, TRANSPORT)

#if IS_SHM_TRANSPORT
#include <omp.h>
#endif

#include "am++/am++.hpp"
#define TRANSPORT_HEADER <am++/BOOST_JOIN(TRANSPORT, _transport).hpp>
#include TRANSPORT_HEADER
#include "am++/basic_coalesced_message_type.hpp"
#include "am++/counter_coalesced_message_type.hpp"
#include "am++/reductions.hpp"
// #include "am++/lock_free_coalesced_message_type.hpp"
// #include "am++/lock_free_coalesced_message_type_packed.hpp"
#include "am++/object_based_addressing.hpp"
#include "am++/message_type_generators.hpp"
#if IS_MPI_TRANSPORT
#include <am++/mpi_sinha_kale_ramkumar_termination_detector.hpp>
#include <mpi.h>
#include <am++/make_mpi_datatype.hpp>
#endif
#include <stdio.h>
#include <string>
// #include <boost/serialization/string.hpp>
// #include <boost/graph/rmat_graph_generator.hpp>
#include "distributed_erdos_renyi_generator.hpp"
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/two_bit_color_map.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/functional/hash.hpp>
#include <boost/assert.hpp>
#include <utility>
#include <functional>
#include <iostream>
#include <sstream>

typedef amplusplus::transport::rank_type rank_type;

typedef boost::compressed_sparse_row_graph<boost::directedS, boost::no_property, boost::no_property, boost::no_property, uint32_t, uint32_t> GraphT;
typedef boost::graph_traits<GraphT>::vertex_descriptor Vertex;

struct update_vertex_data {
  Vertex v;
  size_t dist;

  update_vertex_data(Vertex v, size_t dist): v(v), dist(dist) {}
};

struct get_vertex {
  typedef Vertex result_type;
  Vertex operator()(const update_vertex_data& d) const {return d.v;}
};

inline std::ostream& operator<<(std::ostream& o, const update_vertex_data& d) {
  o << "(" << d.v << " d=" << d.dist << ")";
  return o;
}

struct distrib_type {
  size_t chunk_size;
  explicit distrib_type(size_t chunk_size = 0): chunk_size(chunk_size) {}
  int operator()(size_t x) const {return x / chunk_size;}
};

struct owner_map_type {
  size_t chunk_size;
  explicit owner_map_type(size_t chunk_size = 0): chunk_size(chunk_size) {}
  friend int get(const owner_map_type& o, const update_vertex_data& d) {
    return d.v / o.chunk_size;
  }
};

#if IS_MPI_TRANSPORT
namespace amplusplus {
  template <>
  struct make_mpi_datatype<update_vertex_data> : make_mpi_datatype_base {
    make_mpi_datatype<Vertex> dt1;
    make_mpi_datatype<size_t> dt2;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1(), dt2() {
      int blocklengths[2] = {1, 1};
      MPI_Aint displacements[2];
      update_vertex_data test_object(0, 0);
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.v, &displacements[0]);
      MPI_Get_address(&test_object.dist, &displacements[1]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      MPI_Datatype types[2] = {dt1.get(), dt2.get()};
      MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}
#endif // IS_MPI_TRANSPORT

struct update_vertex_handler {
  Vertex my_start, my_end;
  boost::two_bit_color_map<>* visited;
  std::vector<Vertex>* Q_tail;

  // update_vertex_handler(): my_start(0), my_end(0), visited(NULL), Q_tail(NULL) {}
  update_vertex_handler(Vertex my_start, Vertex my_end, boost::two_bit_color_map<>& visited, std::vector<Vertex>& Q_tail)
    : my_start(my_start), my_end(my_end), visited(&visited), Q_tail(&Q_tail) {}

  // void operator()(rank_type source, const update_vertex_data& data) const {
  void operator()(const update_vertex_data& data) const {
    BOOST_ASSERT (data.v >= my_start && data.v < my_end);
    // fprintf(stderr, "Got %zu %s range\n", (size_t)data.v, (data.v >= my_start && data.v < my_end) ? "in" : "out of");
    if (get(*visited, data.v - my_start) == boost::two_bit_white) {
      // fprintf(stderr, "Enqueueing %zu\n", v);
      put(*visited, data.v - my_start, boost::two_bit_gray);
      Q_tail->push_back(data.v);
    }
#if 0
    // fprintf(stderr, "Got %zu\n", v);
    if (__sync_bool_compare_and_swap(&(*visited)[v - my_start], 0, 1)) {
      // fprintf(stderr, "Enqueueing %zu\n", v);
      *__sync_fetch_and_add(Q_tail, sizeof(size_t)) = v;
    }
#endif
  }
};

template <typename CoalescingTag, typename AMTransport, typename OwnerMap = void> struct coalescing;

struct basic_ct {static const char* print() {return "Basic coalescing";}};
struct counter_ct {static const char* print() {return "Counter coalescing";}};
struct lock_free_ct {static const char* print() {return "Lock-free coalescing";}};
struct lock_free_packed_ct {static const char* print() {return "Lock-free packed coalescing";}};
struct oba_ct {static const char* print() {return "Object-based addressing, direct sends";}};
struct oba_ring_ct {static const char* print() {return "Object-based addressing, ring";}};
struct oba_hypercube_ct {static const char* print() {return "Object-based addressing, hypercube";}};
struct oba_dissemination_ct {static const char* print() {return "Object-based addressing, dissemination";}};
struct oba_hypercube_reduced_ct {static const char* print() {return "Object-based addressing, hypercube, duplicate removal";}};

#if 0
template <typename AMTransport, typename OwnerMap> struct coalescing<basic_ct, AMTransport, OwnerMap> {
  typedef amplusplus::basic_coalesced_message_type<update_vertex_data, update_vertex_handler, AMTransport> type;
  static type make(AMTransport& transport, const OwnerMap&) {return type(transport, (1 << 12));}
};

template <typename AMTransport, typename OwnerMap> struct coalescing<lock_free_packed_ct, AMTransport, OwnerMap> {
  typedef amplusplus::lock_free_coalesced_message_type_packed<update_vertex_data, update_vertex_handler, AMTransport, (1 << 12)> type;
  static type make(AMTransport& transport, const OwnerMap&) {return type(transport);}
};
#endif

template <typename AMTransport, typename OwnerMap>
struct coalescing<oba_ct, AMTransport, OwnerMap> {
  template <typename K>
  static void go(const amplusplus::transport& trans, const OwnerMap& owner, const K& k) {
    typedef amplusplus::simple_generator<amplusplus::basic_coalesced_message_type_gen> Gen;
    Gen gen(amplusplus::basic_coalesced_message_type_gen(1 << 12));
    typename Gen::template call_result<update_vertex_data, update_vertex_handler, OwnerMap, amplusplus::no_reduction_t>::type
      c(gen, trans, owner, amplusplus::no_reduction);
    k(c);
  }
};

template <typename AMTransport, typename OwnerMap> struct coalescing<oba_ring_ct, AMTransport, OwnerMap> {
  template <typename K>
  static void go(AMTransport& transport, const OwnerMap& owner, const K& k) {
    typedef amplusplus::routing_generator<amplusplus::basic_coalesced_message_type_gen, amplusplus::ring_routing> Gen;
    Gen gen(amplusplus::basic_coalesced_message_type_gen(1 << 12), amplusplus::ring_routing(transport.rank(), transport.size()));
    typename Gen::template call_result<update_vertex_data, update_vertex_handler, OwnerMap, amplusplus::no_reduction_t>::type
      c(gen, transport, owner, amplusplus::no_reduction);
    k(c);
  }
};

template <typename AMTransport, typename OwnerMap> struct coalescing<oba_hypercube_ct, AMTransport, OwnerMap> {
  template <typename K>
  static void go(AMTransport& transport, const OwnerMap& owner, const K& k) {
    typedef amplusplus::routing_generator<amplusplus::basic_coalesced_message_type_gen, amplusplus::hypercube_routing> Gen;
    Gen gen(amplusplus::basic_coalesced_message_type_gen(1 << 12), amplusplus::hypercube_routing(transport.rank(), transport.size()));
    typename Gen::template call_result<update_vertex_data, update_vertex_handler, OwnerMap, amplusplus::no_reduction_t>::type
      c(gen, transport, owner, amplusplus::no_reduction);
    k(c);
  }
};

template <typename AMTransport, typename OwnerMap>
struct coalescing<oba_dissemination_ct, AMTransport, OwnerMap> {
  template <typename K>
  static void go(AMTransport& transport, const OwnerMap& owner, const K& k) {
    typedef amplusplus::routing_generator<amplusplus::basic_coalesced_message_type_gen, amplusplus::dissemination_routing> Gen;
    Gen gen(amplusplus::basic_coalesced_message_type_gen(1 << 12), amplusplus::dissemination_routing(transport.rank(), transport.size()));
    typename Gen::template call_result<update_vertex_data, update_vertex_handler, OwnerMap, amplusplus::no_reduction_t>::type
      c(gen, transport, owner, amplusplus::no_reduction);
    k(c);
  }
};

template <typename OwnerMap, typename AMTransport>
struct coalescing<oba_hypercube_reduced_ct, AMTransport, OwnerMap> {
  template <typename K>
  static void go(const amplusplus::transport& trans, const OwnerMap& owner, const K& k) {
    typedef amplusplus::cache_generator<amplusplus::basic_coalesced_message_type_gen, amplusplus::hypercube_routing> Gen;
    Gen gen(amplusplus::basic_coalesced_message_type_gen(1 << 12), 12, amplusplus::hypercube_routing(trans.rank(), trans.size()));
    typename Gen::template call_result<update_vertex_data, update_vertex_handler, OwnerMap, amplusplus::duplicate_removal_t<get_vertex> >::type
      c(gen, trans, owner, amplusplus::duplicate_removal(get_vertex()));
    k(c);
  }
};

struct value_with_loc {long val; int idx;}; // To use MPI::MAXLOC

#if IS_SHM_TRANSPORT
boost::scoped_ptr<boost::barrier> reduction_barrier;
boost::mutex reduction_lock;
value_with_loc reduction_value;
#endif

template <typename Graph, typename AMTransport, typename CT>
struct do_test_2 {
  const Graph& my_graph;
  AMTransport& transport;
  Vertex my_size, my_start, my_end, chunk_size, graph_size;

  do_test_2(const Graph& my_graph, AMTransport& transport, Vertex my_size, Vertex my_start, Vertex my_end, Vertex chunk_size, Vertex graph_size): my_graph(my_graph), transport(transport), my_size(my_size), my_start(my_start), my_end(my_end), chunk_size(chunk_size), graph_size(graph_size) {}

  template <typename UpdateMsg>
  void operator()(UpdateMsg& update_msg) const {
    rank_type rank, size;

    rank = transport.rank();
    size = transport.size();

    std::vector<Vertex> local_queue;
    value_with_loc highest_degree_vertex = {0, 0}; // First is degree, second is vertex index
    BGL_FORALL_VERTICES_T(v, my_graph, Graph) {
      if (out_degree(v, my_graph) > Vertex(highest_degree_vertex.val)) {
        highest_degree_vertex.val = out_degree(v, my_graph);
        highest_degree_vertex.idx = v + my_start;
      }
    }

    boost::two_bit_color_map<> visited(my_size);
    for (Vertex i = 0; i < my_size; ++i) put(visited, i, boost::two_bit_white);

    value_with_loc highest_degree_vertex_global;
#if IS_MPI_TRANSPORT
    MPI_Allreduce(&highest_degree_vertex, &highest_degree_vertex_global, 1, MPI_LONG_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#elif IS_SHM_TRANSPORT
    {
      reduction_barrier->wait();
      if (rank == 0) reduction_value = highest_degree_vertex;
      reduction_barrier->wait();
      {
        boost::lock_guard<boost::mutex> l(reduction_lock);
        if (highest_degree_vertex.val < reduction_value.val) {
          reduction_value = highest_degree_vertex;
        }
      }
      reduction_barrier->wait();
      highest_degree_vertex_global = reduction_value;
      reduction_barrier->wait();
    }
#else
#error "Unhandled transport"
#endif
    if ((Vertex)highest_degree_vertex_global.idx >= my_start &&
        (Vertex)highest_degree_vertex_global.idx < my_end) {
      local_queue.push_back(Vertex(highest_degree_vertex_global.idx));
      put(visited, Vertex(highest_degree_vertex_global.idx) - my_start, boost::two_bit_gray);
    }

    update_msg.set_handler(update_vertex_handler(my_start, my_end, visited, local_queue));

    double start = amplusplus::get_time();
    size_t current_dist = 0;
    size_t local_visited_in_queue = 0;
    while (true) {
      std::vector<Vertex> old_queue;
      old_queue.swap(local_queue);
      std::vector<Vertex>::const_iterator b = old_queue.begin();
      std::vector<Vertex>::const_iterator e = old_queue.end();
      // fprintf(stderr, "%d: On distance %zu, have %zu vertices in queue\n", rank, current_dist, e - b);
      local_visited_in_queue += (e - b);
      // This goes one level too far, but is easier than trying to figure out the correct value
      const unsigned long my_queue_is_empty = old_queue.empty();
      unsigned long all_queues_are_empty;
      {
        // fprintf(stderr, "Begin: Transport = %p\n", &transport);
        amplusplus::scoped_epoch_value epoch(transport, my_queue_is_empty, all_queues_are_empty);
        for (; b != e; ++b) {
          // fprintf(stderr, "Working on vertex %zu\n", *b);
          BGL_FORALL_ADJ_T(Vertex(*b - my_start), w, my_graph, Graph) {
            // fprintf(stderr, "Sending update on %zu to %zu\n", w, w / chunk_size);
#if 0
            if (w >= my_start && w < my_end) {
              if (get(visited, w - my_start) == boost::two_bit_white) {
                put(visited, w - my_start, boost::two_bit_gray);
                local_queue.push_back(w);
              }
            } else {
#endif
              // update_msg.send(update_vertex_data(w, current_dist + 1), w / chunk_size);
              update_msg.send(update_vertex_data(w, current_dist + 1));
#if 0
            }
#endif
          }
        }
        // fprintf(stderr, "End: Transport = %p\n", &transport);
      }
      // if (rank == 0) fprintf(stderr, "all_queues_are_empty = %lu\n", all_queues_are_empty);
      if (all_queues_are_empty == (unsigned long)size) break;
      ++current_dist;
    }
    double stop = amplusplus::get_time();

    unsigned long local_visited = 0;
    for (Vertex i = 0; i < my_size; ++i) {
      if (get(visited, i) != boost::two_bit_white) ++local_visited;
    }

    unsigned long global_visited = 0;
    unsigned long global_visited_in_queue = 0;
    {
      amplusplus::scoped_epoch_value epoch(transport, local_visited, global_visited);
    }
    {
      amplusplus::scoped_epoch_value epoch(transport, local_visited_in_queue, global_visited_in_queue);
    }
    if (rank == 0) fprintf(stderr, "%s: %zu vertices (%zu from queue) visited of %zu\n", CT::print(), (size_t)global_visited, (size_t)global_visited_in_queue, (size_t)graph_size);
    if (rank == 0) fprintf(stdout, "%s: BFS on %zu vertices took %lf s on %zu procs (%zu dist)\n", CT::print(), (size_t)graph_size, stop - start, size, current_dist);
  }
};

template <typename Graph, typename AMTransport, typename CT>
void do_test(const Graph& my_graph, AMTransport& transport, CT, Vertex my_size, Vertex my_start, Vertex my_end, Vertex chunk_size, Vertex graph_size) {
  coalescing<CT, AMTransport, owner_map_type>::go(transport, owner_map_type(chunk_size), do_test_2<Graph, AMTransport, CT>(my_graph, transport, my_size, my_start, my_end, chunk_size, graph_size));
}

void do_one_thread(amplusplus::environment& env) {
  // fprintf(stderr, "do_one_thread on %u\n", tid);

  amplusplus::transport trans = env.create_transport();
  rank_type rank = trans.rank();
  rank_type size = trans.size();
  // fprintf(stderr, "have transport in %u\n", tid);

#if IS_MPI_TRANSPORT
  amplusplus::register_mpi_datatype<update_vertex_data>();
#endif

  const Vertex graph_size = 1000000 * size;
  const Vertex chunk_size = (graph_size + size - 1) / size;
  const Vertex my_start = chunk_size * rank;
  const Vertex my_end = (std::min)(chunk_size * (rank + 1), rank_type(graph_size));
  const Vertex my_size = my_end - my_start;

  std::vector<std::pair<uint32_t, uint32_t> > edges;

  const size_t num_edges = 8 * graph_size;
  if (rank == 0) fprintf(stderr, "Testing with %zu vertices and %zu edges\n", (size_t)graph_size, (size_t)num_edges);

  // typedef boost::rand48 Generator;
  typedef boost::random::splittable_ecuyer1988 Generator;
  Generator gen;
  // boost::rmat_iterator<Generator, GraphT> edges_b(gen, graph_size, num_edges, .54, .21, .21, .04, true);
  // boost::rmat_iterator<Generator, GraphT> edges_b(gen, graph_size, num_edges, 0.57, 0.19, 0.19, 0.05, true);
  // boost::rmat_iterator<Generator, GraphT> edges_b(gen, graph_size, num_edges, 0.61, 0.17, 0.17, 0.05, true);
  // boost::rmat_iterator<Generator, GraphT> edges_e;
  boost::distributed_erdos_renyi_iterator<distrib_type, Generator, Vertex>
    edges_b(gen, graph_size, double(num_edges) / graph_size / graph_size, false, false, distrib_type(chunk_size), rank),
    edges_e;
  if (rank == 0) fprintf(stderr, "Starting to build edge lists\n");

  unsigned long edge_count = 0, total_edge_count = 0;
  for (; edges_b != edges_e; ++edges_b) {
    std::pair<Vertex, Vertex> e = *edges_b;
    // fprintf(stderr, "Edge %lu %lu of %lu-%lu\n", (unsigned long)e.first, (unsigned long)e.second, (unsigned long)my_start, (unsigned long)my_end);
    if (e.first >= my_start && e.first < my_end && e.second < graph_size) {
      edges.push_back(std::make_pair(e.first - my_start, e.second));
    }
    ++edge_count;
    if (rank == 0) {
      if (edge_count % 10000000 == 0) fprintf(stderr, "On edge %lu of %zu\n", (unsigned long)edge_count, (size_t)(num_edges / size));
    }
  }
  {
    amplusplus::scoped_epoch_value epoch(trans, edge_count, total_edge_count);
  }
  if (rank == 0) fprintf(stderr, "Made %lu edges total\n", (unsigned long)edge_count);
  // fprintf(stderr, "Adding Hamiltonian cycle\n");
  boost::variate_generator<Generator, boost::uniform_int<> > random_vertex(gen, boost::uniform_int<>(0, graph_size - 1));
  std::vector<Vertex> perm(graph_size);
  for (Vertex i = 0; i < graph_size; ++i) perm[i] = i;
  std::random_shuffle(perm.begin(), perm.end(), random_vertex);
  size_t raw_edge_count = edges.size();
  for (Vertex i = 0; i < graph_size; ++i) {
    if (perm[i] >= my_start && perm[i] < my_end) {
      edges.push_back(std::make_pair(perm[i] - my_start, perm[(i + 1) % graph_size]));
    }
  }
  std::sort(edges.begin() + raw_edge_count, edges.end());
  std::inplace_merge(edges.begin(), edges.begin() + raw_edge_count, edges.end());
  // if (rank == 0) fprintf(stderr, "Sorting edge list\n");
  // std::sort(edges.begin(), edges.end());
  if (rank == 0) fprintf(stderr, "Building graph\n");

  // This graph has out edges pointing to nonexistent vertices, but this works
  // for what we need.  The source vertex indexes are biased by -my_start.
  GraphT my_graph(boost::edges_are_sorted, edges.begin(), edges.end(), my_size);
  if (rank == 0) fprintf(stderr, "Running BFS\n");

#if 0
  {
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(make_mpi_sinha_kale_ramkumar_termination_detector(transport));
    do_test(my_graph, transport, basic_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }
#endif

  {
    if (rank == 0) std::cout << "SKR termination detector (test only)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    // transport.set_termination_detector(make_mpi_sinha_kale_ramkumar_termination_detector(transport));
    do_test(my_graph, transport, oba_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

  {
    if (rank == 0) std::cout << "SKR termination detector (no routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    // transport.set_termination_detector(make_mpi_sinha_kale_ramkumar_termination_detector(transport));
    do_test(my_graph, transport, oba_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

#if 0
#if IS_MPI_TRANSPORT
  {
    if (rank == 0) std::cout << "PCX termination detector (no routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_pcx_termination_detector(transport));
    do_test(my_graph, transport, oba_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

  {
    if (rank == 0) std::cout << "NBX termination detector (no routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_nbx_termination_detector(transport));
    do_test(my_graph, transport, oba_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }
#endif
#endif

#if 0
  {
    amplusplus::transport transport = env.create_transport();
    do_test(my_graph, transport, oba_ring_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }
#endif

  if ((size & (size - 1)) == 0) { // Is power of 2
    if (rank == 0) std::cout << "SKR termination detector (hypercube routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    // transport.set_termination_detector(make_mpi_sinha_kale_ramkumar_termination_detector(transport));
    do_test(my_graph, transport, oba_hypercube_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

#if 0
#if IS_MPI_TRANSPORT
  if ((size & (size - 1)) == 0) { // Is power of 2
    if (rank == 0) std::cout << "NBX termination detector (hypercube routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_nbx_termination_detector(transport));
    do_test(my_graph, transport, oba_hypercube_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

  if ((size & (size - 1)) == 0) { // Is power of 2
    if (rank == 0) std::cout << "NBX termination detector (hypercube routing, reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_nbx_termination_detector(transport));
    do_test(my_graph, transport, oba_hypercube_reduced_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

  if ((size & (size - 1)) == 0) { // Is power of 2
    if (rank == 0) std::cout << "PCX termination detector (hypercube routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_pcx_termination_detector(transport));
    do_test(my_graph, transport, oba_hypercube_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

  if ((size & (size - 1)) == 0) { // Is power of 2
    if (rank == 0) std::cout << "PCX termination detector (hypercube routing, reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_pcx_termination_detector(transport));
    do_test(my_graph, transport, oba_hypercube_reduced_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }
#endif

  {
    if (rank == 0) std::cout << "SKR termination detector (dissemination routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    do_test(my_graph, transport, oba_dissemination_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }

#if IS_MPI_TRANSPORT
  {
    if (rank == 0) std::cout << "NBX termination detector (dissemination routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_nbx_termination_detector(transport));
    do_test(my_graph, transport, oba_dissemination_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }
  {
    if (rank == 0) std::cout << "PCX termination detector (dissemination routing, no reduction)" << std::endl;
    amplusplus::transport transport = env.create_transport();
    transport.set_termination_detector(amplusplus::make_mpi_pcx_termination_detector(transport));
    do_test(my_graph, transport, oba_dissemination_ct(), my_size, my_start, my_end, chunk_size, graph_size);
  }
#endif
#endif

}

int main(int argc, char** argv) {
  (void)argc; (void)argv;
#if IS_MPI_TRANSPORT
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  do_one_thread(env);
#elif IS_SHM_TRANSPORT
  boost::scoped_ptr<amplusplus::shm_environment_common> common;

#pragma omp parallel
  {
#pragma omp single
    {
      common.reset(new amplusplus::shm_environment_common(omp_get_num_threads()));
      reduction_barrier.reset(new boost::barrier(omp_get_num_threads()));
    }
    amplusplus::environment env = amplusplus::shm_environment(*common, omp_get_thread_num());
    do_one_thread(env);
  }
#endif
  return 0;
}
