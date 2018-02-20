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
#include "am++/counter_coalesced_message_type.hpp"
#include "am++/reductions.hpp"
#include "am++/message_type_generators.hpp"
#if IS_MPI_TRANSPORT
#include <mpi.h>
#endif
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <stdio.h>
#include <string>
#include <boost/serialization/string.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <utility>
#include <functional>

typedef amplusplus::transport::rank_type rank_type;

typedef std::pair<size_t, double> accum_b_data;

struct accum_b_handler {
  std::vector<double>* b;
  size_t my_start;
  size_t* count;
  accum_b_handler(): b(NULL), my_start(0), count(NULL) {}
  accum_b_handler(std::vector<double>& b, size_t my_start, size_t& count): b(&b), my_start(my_start), count(&count) {}

  void operator()(const accum_b_data& data) const {
    // fprintf(stderr, "%d got %zu -> %d from %d\n", rank, i, val, source);
    // __sync_fetch_and_add(&b[data.first - my_start], data.second);
    // __sync_fetch_and_add(&count, 1);
    (*b)[data.first - my_start] += data.second;
    ++(*count);
  }
};

struct double_sum_binop {
  typedef size_t key_type;
  typedef double value_type;
  typedef std::pair<key_type, value_type> pair_type;
  key_type& get_key(pair_type& p) const {return p.first;}
  const key_type& get_key(const pair_type& p) const {return p.first;}
  value_type& get_value(pair_type& p) const {return p.second;}
  const value_type& get_value(const pair_type& p) const {return p.second;}
  pair_type make_pair(const key_type& k, const value_type& v) const {return std::make_pair(k, v);}
  size_t hash_key(key_type k) const {return k;}
  key_type dummy_key(size_t h) const {return h + 1;}
  bool is_valid_key(size_t h, key_type k) const {return k != h + 1;}
  double op(double a, double b) const {return a + b;}
  bool is_identity(double x) const {return x == 0;}
};

struct owner_map_type {
  size_t chunk_size;
  explicit owner_map_type(size_t chunk_size = 0): chunk_size(chunk_size) {}
  friend int get(const owner_map_type& o, const accum_b_data& d) {
    return d.first / o.chunk_size;
  }
};

template <typename Graph>
void do_test(amplusplus::transport trans, const Graph& my_graph, const std::vector<double>& x_local, std::vector<double>& b_local, const std::vector<double>& A_local, size_t /*my_size*/, size_t my_start, size_t /*my_end*/, size_t chunk_size, size_t graph_size) {
  rank_type rank = trans.rank();
  rank_type size = trans.size();

  size_t count = 0;

#if IS_MPI_TRANSPORT
  amplusplus::register_mpi_datatype<accum_b_data>();
#endif

  BGL_FORALL_VERTICES_T(v, my_graph, Graph) {
    b_local[v] = 0;
  }

  {amplusplus::scoped_epoch e(trans);}
  double start = amplusplus::get_time();
  {
  typedef amplusplus::cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::hypercube_routing> Gen;
  Gen gen(amplusplus::counter_coalesced_message_type_gen(1 << 12), 9, amplusplus::hypercube_routing(trans.rank(), trans.size()));
  owner_map_type owner(chunk_size);
  typedef typename Gen::template call_result<accum_b_data, accum_b_handler, owner_map_type, amplusplus::combination_t<std::plus<double>, double> >::type accum_b_type;
  accum_b_type
    accum_b_msg(gen, trans, owner, amplusplus::combination(std::plus<double>(), 0.));
  accum_b_msg.set_handler(accum_b_handler(b_local, my_start, count));

    amplusplus::scoped_epoch epoch(trans);
    BGL_FORALL_EDGES_T(e, my_graph, Graph) {
      size_t idx = get(boost::edge_index, my_graph, e);
      size_t t = target(e, my_graph);
      double coeff = A_local[idx];
      accum_b_msg.send(accum_b_data(t, coeff * x_local[source(e, my_graph)]));
      // if (idx % 10000 == 0) fprintf(stderr, "On edge %zu\n", (size_t)idx);
    }
  }
  double stop = amplusplus::get_time();

  fprintf(stderr, "Rank %zu did %zu additions and has %zu edges.\n", rank, count, num_edges(my_graph));
  if (rank == 0) fprintf(stdout, "Matvec_t on %zu vertices took %lf s on %zu procs\n", graph_size, stop - start, size);
}

void do_one_thread(amplusplus::environment& env) {
  rank_type rank, size;
  {
    amplusplus::transport trans = env.create_transport();
    rank = trans.rank();
    size = trans.size();
  }

  const size_t graph_size = 40000; //  * size;
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

  for (size_t i = 0; i < my_size; ++i) source_array[i] = (double)(drand48() * 100);
  for (size_t i = 0; i < boost::num_edges(my_graph); ++i) coeffs[i] = (double)(drand48() * 100);

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
