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
#include <boost/config.hpp>
#define TRANSPORT_HEADER <am++/BOOST_JOIN(TRANSPORT, _transport).hpp>
#include TRANSPORT_HEADER
#include "am++/basic_coalesced_message_type.hpp"
#if IS_MPI_TRANSPORT
#include <mpi.h>
#endif
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <stdio.h>
#include <string>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/assert.hpp>
#include <utility>
#include <functional>

typedef amplusplus::transport::rank_type rank_type;

template <typename K, typename V, typename AMTransport>
struct read_only_am_property_map : boost::noncopyable {
  struct get_handler {
    read_only_am_property_map* pm;
    get_handler(): pm(NULL) {}
    get_handler(read_only_am_property_map& pm): pm(&pm) {}
    void operator()(rank_type src, const K& key) const {
      // fprintf(stderr, "%zu: Get request for %zu from %zu\n", pm->rank, key, src);
      pm->put_ghost_msg.send(std::make_pair(key, pm->local_data[key - pm->my_start]), src);
    }
  };

  struct put_ghost_handler {
    read_only_am_property_map* pm;
    put_ghost_handler(): pm(NULL) {}
    put_ghost_handler(read_only_am_property_map& pm): pm(&pm) {}
    void operator()(rank_type src, const std::pair<K, V>& req) const {
      (void)src;
      // fprintf(stderr, "%zu: Put_ghost request for %zu from %zu\n", pm->rank, req.first, src);
      pm->ghost_cells[req.first] = req.second;
    }
  };

  const int dummy_first_member_for_init_order; // Unused
  rank_type rank;
  const std::vector<V>& local_data;
  K my_start;
  K chunk_size;
  std::map<K, V> ghost_cells;
  AMTransport& transport;
  typedef amplusplus::basic_coalesced_message_type<K, get_handler> get_type;
  typedef amplusplus::basic_coalesced_message_type<std::pair<K, V>, put_ghost_handler> put_ghost_type;
  get_type get_msg;
  put_ghost_type put_ghost_msg;

  read_only_am_property_map(AMTransport& transport, const std::vector<V>& local_data, K my_start, K chunk_size):
    dummy_first_member_for_init_order((
#if IS_MPI_TRANSPORT
      amplusplus::register_mpi_datatype<std::pair<K, V> >(),
#endif
      0)),
    rank(transport.rank()), local_data(local_data), my_start(my_start), chunk_size(chunk_size),
    ghost_cells(), transport(transport), get_msg(amplusplus::basic_coalesced_message_type_gen(1 << 10), transport), put_ghost_msg(amplusplus::basic_coalesced_message_type_gen(1 << 10), transport)
  {
    // fprintf(stderr, "Registering handlers\n");
    get_msg.set_handler(get_handler(*this));
    put_ghost_msg.set_handler(put_ghost_handler(*this));
  }

  V get_local_val(K key) const {return local_data[key - my_start];}
  void start_remote_get(K key) {
    if (ghost_cells.find(key) != ghost_cells.end()) return;
    // ghost_cells[key] = 0. / 0.; // Set to NaN so it is marked as already sent
    // fprintf(stderr, "%zu: Sending get to %d for %zu\n", rank, (int)(key / chunk_size), key);
    get_msg.send(key, key / chunk_size);
  }
  V get_remote_val(K key) const {
#if 0
    if (ghost_cells.find(key) == ghost_cells.end()) {
      fprintf(stderr, "%zu: Failed to find %zu sent by %d\n", rank, key, (int)(key / chunk_size));
    }
#endif
    BOOST_ASSERT (ghost_cells.find(key) != ghost_cells.end());
    return ghost_cells.find(key)->second;
  }
};

template <typename Graph, typename AMTransport>
void do_test(AMTransport& transport, const Graph& my_graph, const std::vector<double>& x_local, std::vector<double>& b_local, const std::vector<double>& A_local, size_t my_size, size_t my_start, size_t my_end, size_t chunk_size, size_t graph_size) {
  rank_type rank = transport.rank();
  rank_type size = transport.size();

  typedef read_only_am_property_map<size_t, double, AMTransport> pm_type;
  pm_type pm(transport, x_local, my_start, chunk_size);
  {amplusplus::scoped_epoch e(transport);}

  double start = amplusplus::get_time();
  {
    amplusplus::scoped_epoch epoch(transport);
    BGL_FORALL_EDGES_T(e, my_graph, Graph) {
      size_t t = target(e, my_graph);
      if (t < my_start || t >= my_end) {
        pm.start_remote_get(t);
      }
    }
  }
  BGL_FORALL_VERTICES_T(u, my_graph, Graph) {
    if (u >= my_size) continue;
    double sum = 0;
    BGL_FORALL_OUTEDGES_T(u, e, my_graph, Graph) {
      size_t v = target(e, my_graph);
      // BOOST_ASSERT (get(boost::edge_index, my_graph, e) < A_local.size());
      double coeff = A_local[get(boost::edge_index, my_graph, e)];
      if (v < my_start || v >= my_end) {
        sum += coeff * pm.get_remote_val(v);
      } else {
        // BOOST_ASSERT (v >= my_start && v < my_end);
        sum += coeff * pm.get_local_val(v);
      }
    }
    // BOOST_ASSERT (u < my_size);
    b_local[u] = sum;
  }
  {amplusplus::scoped_epoch e(transport);}
  double stop = amplusplus::get_time();
  fprintf(stderr, "Rank %zu has %zu ghost cell(s)\n", rank, pm.ghost_cells.size());
  if (rank == 0) fprintf(stdout, "Matvec on %zu vertices took %lf s on %zu procs\n", graph_size, stop - start, size);
}

void do_one_thread(amplusplus::environment& env) {
  rank_type rank, size;
  {
    amplusplus::transport trans = env.create_transport();
    rank = trans.rank();
    size = trans.size();
  }

  const size_t graph_size = 10000 * size;
  const size_t chunk_size = (graph_size + size - 1) / size;
  const size_t my_start = chunk_size * rank;
  const size_t my_end = (std::min)(chunk_size * (rank + 1), graph_size);
  const size_t my_size = my_end - my_start;

  std::vector<std::pair<size_t, size_t> > edges;

  const size_t num_edges = 8 * graph_size;
  if (rank == 0) fprintf(stderr, "Testing with %zu vertices and %zu edges\n", graph_size, num_edges);

  typedef boost::compressed_sparse_row_graph<> Graph;
  typedef boost::rand48 Generator;
  Generator gen;
  boost::rmat_iterator<Generator, Graph> edges_b(gen, graph_size, num_edges, .54, .21, .21, .04, true);
  boost::rmat_iterator<Generator, Graph> edges_e;
  if (rank == 0) fprintf(stderr, "Starting to build edge lists\n");

  unsigned int edge_count = 0;
  for (; edges_b != edges_e; ++edges_b) {
    if (edges_b->first >= my_start && edges_b->first < my_end) {
      edges.push_back(std::make_pair(edges_b->first - my_start, edges_b->second));
    }
    if (rank == 0) {
      ++edge_count;
      if (edge_count % 1000000 == 0) fprintf(stderr, "On edge %u of %zu\n", edge_count, num_edges);
    }
  }
  if (rank == 0) fprintf(stderr, "Sorting edge list\n");
  std::sort(edges.begin(), edges.end());
  if (rank == 0) fprintf(stderr, "Building graph\n");

  // This graph has out edges pointing to nonexistent vertices, but this works
  // for what we need.  The source vertex indexes are biased by -my_start.
  boost::compressed_sparse_row_graph<> my_graph(boost::edges_are_sorted, edges.begin(), edges.end(), my_size);

  if (rank == 0) fprintf(stderr, "Initializing coeffs\n");
  std::vector<double> source_array(my_size), dest_array(my_size);
  std::vector<double> coeffs(boost::num_edges(my_graph));

  for (size_t i = 0; i < my_size; ++i) source_array[i] = drand48();
  for (size_t i = 0; i < boost::num_edges(my_graph); ++i) coeffs[i] = drand48();

  if (rank == 0) fprintf(stderr, "Running matvec\n");

  {
    amplusplus::transport transport = env.create_transport();
    do_test(transport, my_graph, source_array, dest_array, coeffs, my_size, my_start, my_end, chunk_size, graph_size);
  }
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
    }
    amplusplus::environment env = amplusplus::shm_environment(*common, omp_get_thread_num());
    do_one_thread(env);
  }
#endif
  return 0;
}
