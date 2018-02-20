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

#include <am++/am++.hpp>
// #include <am++/lock_free_coalesced_message_type.hpp>
// #include <am++/lock_free_coalesced_message_type_packed.hpp>
#include <am++/basic_coalesced_message_type.hpp>
#include <am++/counter_coalesced_message_type.hpp>
#include <am++/reductions.hpp>
#include <am++/message_type_generators.hpp>
#include <am++/detail/thread_support.hpp>
#include <boost/config.hpp>
#define TRANSPORT_HEADER <am++/BOOST_JOIN(TRANSPORT, _transport).hpp>
#include TRANSPORT_HEADER
#if IS_MPI_TRANSPORT
#include <mpi.h>
#include <am++/make_mpi_datatype.hpp>
#endif
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/type_traits.hpp>
#include <stdio.h>
#include <string>
// #include <boost/serialization/string.hpp>
// #include <boost/graph/rmat_graph_generator.hpp>
#include "rmat_graph_generator_faster.hpp"
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/two_bit_color_map.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <utility>
#include <functional>

typedef boost::compressed_sparse_row_graph<boost::directedS, boost::no_property, boost::no_property, boost::no_property, uint32_t, uint32_t> GraphT;
typedef boost::graph_traits<GraphT>::vertex_descriptor Vertex;
typedef amplusplus::rank_type rank_type;

static const int buffer_size = (1 << 10);

struct update_vertex_data {
  Vertex v;
  size_t dist;

  update_vertex_data(): v(0), dist(0) {}
  update_vertex_data(Vertex v, size_t dist): v(v), dist(dist) {}
};

struct get_vertex {
  typedef Vertex result_type;
  Vertex operator()(const update_vertex_data& d) const {return d.v;}
};

struct update_vertex_data_remove_duplicate_policy {
  typedef update_vertex_data value_type;
  typedef Vertex stored_type;
  Vertex stored_value(const update_vertex_data& d) const {
    return d.v;
  }
  bool match(const update_vertex_data& d, Vertex stored) const {
    return d.v == stored;
  }
  size_t hash(Vertex v) const {
    return (size_t)v;
  }
  Vertex dummy_value(int h) const {
    return Vertex(h + 1);
  }
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
#endif

struct update_vertex_handler {
  Vertex my_start, my_end;
  char* visited;
  Vertex** local_queue;
  amplusplus::detail::atomic<size_t>* local_queue_size;

  update_vertex_handler(): my_start(0), my_end(0), visited(NULL), local_queue(NULL), local_queue_size(NULL) {}
  update_vertex_handler(Vertex my_start, Vertex my_end, char* visited, Vertex*& local_queue, amplusplus::detail::atomic<size_t>& local_queue_size)
    : my_start(my_start), my_end(my_end), visited(visited), local_queue(&local_queue), local_queue_size(&local_queue_size) {}
#ifndef BOOST_NO_DEFAULTED_FUNCTIONS
  update_vertex_handler(update_vertex_handler&&) = default;
  update_vertex_handler(const update_vertex_handler&) = delete;
  update_vertex_handler& operator=(update_vertex_handler&&) = default;
  update_vertex_handler& operator=(const update_vertex_handler&) = delete;
#endif

  void operator()(const update_vertex_data& data) const {
    // fprintf(stderr, "Got %zu %s range\n", v, (v >= my_start && v < my_end) ? "in" : "out of");
    if (__sync_val_compare_and_swap(&visited[data.v - my_start], 0, 1) == 0) {
      // fprintf(stderr, "Enqueueing %zu\n", v);
      (*local_queue)[(*local_queue_size).fetch_add(1)] = data.v;
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

template <typename Graph, typename UpdateMsgType>
struct thread_body {
  boost::barrier& bar;
  amplusplus::transport trans;
  Vertex*& local_queue;
  amplusplus::detail::atomic<size_t>& local_queue_size;
  Vertex*& old_queue;
  size_t& old_queue_size;
  const Graph& my_graph;
  size_t& current_dist;
  size_t& local_visited_in_queue;
  const unsigned int num_threads;
  char* const visited;
  const Vertex my_start;
  const Vertex my_end;
  const Vertex chunk_size;
  UpdateMsgType& update_msg;

  thread_body(boost::barrier& bar, amplusplus::transport trans, Vertex*& local_queue, amplusplus::detail::atomic<size_t>& local_queue_size, Vertex*& old_queue, size_t& old_queue_size, const Graph& my_graph, size_t& current_dist, size_t& local_visited_in_queue, unsigned int num_threads, char* const visited, const Vertex my_start, const Vertex my_end, const Vertex chunk_size, UpdateMsgType& update_msg)
    : bar(bar), trans(trans), local_queue(local_queue), local_queue_size(local_queue_size), old_queue(old_queue), old_queue_size(old_queue_size), my_graph(my_graph), current_dist(current_dist), local_visited_in_queue(local_visited_in_queue), num_threads(num_threads), visited(visited), my_start(my_start), my_end(my_end), chunk_size(chunk_size), update_msg(update_msg)
    {}

  void operator()(int thread_num) {
    rank_type size = trans.size();
    AMPLUSPLUS_WITH_THREAD_ID(thread_num) {
      while (true) {
        bar.wait();
        if (thread_num == 0) {
          std::swap(old_queue, local_queue);
          old_queue_size = local_queue_size.load();
          local_queue_size.store(0);
        }
        bar.wait();
        const Vertex* b = old_queue;
        size_t queue_size = old_queue_size;
        // fprintf(stderr, "%d: On distance %zu, have %zu vertices in queue\n", rank, current_dist, queue_size);
        if (thread_num == 0) {
          local_visited_in_queue += queue_size;
        }
        // This goes one level too far, but is easier than trying to figure out the correct value
        const unsigned long my_queue_is_empty = (thread_num == 0 && queue_size == 0);
        unsigned long all_queues_are_empty;
        {
          amplusplus::scoped_epoch_value epoch(trans, my_queue_is_empty, all_queues_are_empty);
          for (size_t i = thread_num; i < queue_size; i += num_threads) {
            // fprintf(stderr, "Working on vertex %zu\n", *b);
            BGL_FORALL_ADJ_T(Vertex(b[i] - my_start), w, my_graph, Graph) {
              // fprintf(stderr, "Sending update on %zu to %zu\n", w, w / chunk_size);
#if 0
              if (w >= my_start && w < my_end) {
                if (__sync_val_compare_and_swap(&visited[w - my_start], 0, 1) == 0) {
                  local_queue[local_queue_size.fetch_add(1)] = w;
                }
              } else {
#endif
                update_msg.send_with_tid(update_vertex_data(w, current_dist + 1), thread_num);
#if 0
              }
#endif
            }
          }
        }
        // fprintf(stderr, "Thread %ld returning %d to epoch\n", thread_num, (thread_num == 0 && queue_size == 0));
        // if (thread_num == 0 && trans.rank() == 0) fprintf(stderr, "%d Thread %ld: queue_size = %lu, all_queues_are_empty = %lu\n", (int)current_dist, thread_num, queue_size, all_queues_are_empty);
        if (all_queues_are_empty == size) break;
        bar.wait();
        if (thread_num == 0) {
          ++current_dist;
        }
      }
    }
  }
};

struct reduction_none {};
template <size_t Size> struct reduction_simple_cache {};
template <size_t Size> struct reduction_per_thread_cache {};
struct reduction_unordered_set {};

struct value_with_loc {long val; int idx;}; // To use MPI::MAXLOC

#if IS_SHM_TRANSPORT
boost::scoped_ptr<boost::barrier> reduction_barrier;
boost::mutex reduction_lock;
value_with_loc reduction_value;
#endif

template <typename Graph, typename Gen>
void do_test(const Graph& my_graph, amplusplus::transport trans, Vertex my_size, Vertex my_start, Vertex my_end, Vertex chunk_size, Vertex graph_size, int nthreads, const Gen& gen) {
  boost::scoped_array<Vertex> local_queue_storage(new Vertex[my_size]);
  Vertex* local_queue = local_queue_storage.get();
  amplusplus::detail::atomic<size_t> local_queue_size(0);

  rank_type rank = trans.rank();
  rank_type size = trans.size();

  value_with_loc highest_degree_vertex = {0, 0}; // First is degree, second is vertex index
  BGL_FORALL_VERTICES_T(v, my_graph, Graph) {
    if (out_degree(v, my_graph) > (unsigned int)highest_degree_vertex.val) {
      highest_degree_vertex.val = out_degree(v, my_graph);
      highest_degree_vertex.idx = v + my_start;
    }
  }

  boost::scoped_array<char> visited(new char[my_size]);
  for (Vertex i = 0; i < my_size; ++i) visited[i] = 0;

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
    local_queue[local_queue_size.fetch_add(1)] = highest_degree_vertex_global.idx;
    visited[highest_degree_vertex_global.idx - my_start] = 1;
  }

  unsigned int num_threads = nthreads;
  boost::barrier bar(num_threads);
  trans.set_nthreads(num_threads);

  owner_map_type owner(chunk_size);
  typedef typename Gen::template call_result<update_vertex_data, update_vertex_handler, owner_map_type, amplusplus::duplicate_removal_t<get_vertex> >::type update_msg_type;
  update_msg_type update_msg(gen, trans, owner, amplusplus::duplicate_removal(get_vertex()));
  update_msg.set_handler(update_vertex_handler(my_start, my_end, visited.get(), local_queue, local_queue_size));

  double start = amplusplus::get_time();
  size_t current_dist = 0;
  size_t local_visited_in_queue = 0;
  boost::scoped_array<Vertex> old_queue_storage(new Vertex[my_size]);
  Vertex* old_queue = old_queue_storage.get();
  size_t old_queue_size = 0;

  trans.set_nthreads(num_threads);
  boost::scoped_array<boost::thread> threads(new boost::thread[num_threads - 1]);
  for (unsigned int i = 0; i + 1 < num_threads; ++i) {
    boost::thread thr(thread_body<Graph, typename boost::remove_reference<update_msg_type>::type>(bar, trans, local_queue, local_queue_size, old_queue, old_queue_size, my_graph, current_dist, local_visited_in_queue, num_threads, visited.get(), my_start, my_end, chunk_size, update_msg), i);
    threads[i].swap(thr);
  }
  thread_body<Graph, typename boost::remove_reference<update_msg_type>::type>(bar, trans, local_queue, local_queue_size, old_queue, old_queue_size, my_graph, current_dist, local_visited_in_queue, num_threads, visited.get(), my_start, my_end, chunk_size, update_msg)(num_threads - 1);
  for (unsigned int i = 0; i + 1 < num_threads; ++i) {
    threads[i].join();
  }
  double stop = amplusplus::get_time();
  trans.set_nthreads(1);

  unsigned long local_visited = 0;
  for (Vertex i = 0; i < my_size; ++i) {
    if (visited[i] != 0) ++local_visited;
  }

  unsigned long global_visited = 0;
  unsigned long global_visited_in_queue = 0;
  AMPLUSPLUS_WITH_THREAD_ID(0) {
    {
      amplusplus::scoped_epoch_value epoch(trans, local_visited, global_visited);
    }
    {
      amplusplus::scoped_epoch_value epoch(trans, local_visited_in_queue, global_visited_in_queue);
    }
  }
  if (rank == 0) fprintf(stderr, "%lu vertices (%lu from queue) visited of %zu\n", global_visited, global_visited_in_queue, (size_t)graph_size);
  if (rank == 0) fprintf(stdout, "BFS on %zu vertices took %lf s on %zu procs (%zu dist)\n", (size_t)graph_size, stop - start, size, current_dist);
}

int do_one_thread(int argc, char** argv, amplusplus::environment& env) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " nthreads" << std::endl;
    return 1;
  }

  int num_threads = boost::lexical_cast<int>(argv[1]);

  rank_type rank, size;
  {
    amplusplus::transport trans = env.create_transport();
    rank = trans.rank();
    size = trans.size();
  }

#if IS_MPI_TRANSPORT
  amplusplus::register_mpi_datatype<update_vertex_data>();
#endif

  const Vertex graph_size = 100000 * size;
  // const Vertex graph_size = 50000000;
  // const Vertex graph_size = 5000 * size;
  const Vertex chunk_size = (graph_size + size - 1) / size;
  const Vertex my_start = chunk_size * rank;
  const Vertex my_end = (std::min)(chunk_size * Vertex(rank + 1), graph_size);
  const Vertex my_size = my_end - my_start;

  std::vector<std::pair<uint32_t, uint32_t> > edges;

  const size_t num_edges = 8 * graph_size;
  if (rank == 0) fprintf(stderr, "Testing with %zu vertices and %zu edges\n", (size_t)graph_size, (size_t)num_edges);

  typedef boost::rand48 Generator;
  Generator gen;
  // boost::rmat_iterator<Generator, GraphT> edges_b(gen, graph_size, num_edges, .54, .21, .21, .04, true);
  boost::rmat_iterator_faster<Generator, GraphT> edges_b(gen, graph_size, num_edges, 0.57, 0.19, 0.19, 0.05, true);
  // boost::rmat_iterator<Generator, GraphT> edges_b(gen, graph_size, num_edges, 0.61, 0.17, 0.17, 0.05, true);
  boost::rmat_iterator_faster<Generator, GraphT> edges_e;
  if (rank == 0) fprintf(stderr, "Starting to build edge lists\n");

  unsigned int edge_count = 0;
  for (; edges_b != edges_e; ++edges_b) {
    if (edges_b->first >= my_start && edges_b->first < my_end) {
      edges.push_back(std::make_pair(edges_b->first - my_start, edges_b->second));
    }
    if (rank == 0) {
      ++edge_count;
      if (edge_count % 10000000 == 0) fprintf(stderr, "On edge %u of %zu\n", edge_count, num_edges);
    }
  }
  fprintf(stderr, "Adding Hamiltonian cycle\n");
  boost::variate_generator<Generator, boost::uniform_int<> > random_vertex(gen, boost::uniform_int<>(0, graph_size - 1));
  std::vector<Vertex> perm(graph_size);
  for (Vertex i = 0; i < graph_size; ++i) perm[i] = i;
  std::random_shuffle(perm.begin(), perm.end(), random_vertex);
  for (Vertex i = 0; i < graph_size; ++i) {
    if (perm[i] >= my_start && perm[i] < my_end) {
      edges.push_back(std::make_pair(perm[i] - my_start, perm[(i + 1) % graph_size]));
    }
  }
  if (rank == 0) fprintf(stderr, "Sorting edge list\n");
  std::sort(edges.begin(), edges.end());
  if (rank == 0) fprintf(stderr, "Building graph\n");

  // This graph has out edges pointing to nonexistent vertices, but this works
  // for what we need.  The source vertex indexes are biased by -my_start.
  GraphT my_graph(boost::edges_are_sorted, edges.begin(), edges.end(), my_size);
  if (rank == 0) fprintf(stderr, "Running BFS\n");

  {
    if (rank == 0) fprintf(stderr, "Shared transport, no duplicate removal\n");
    amplusplus::transport trans = env.create_transport();
    for (int i = 0; i < 10; ++i) {
      do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen>(amplusplus::counter_coalesced_message_type_gen(1 << 12)));
    }
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, no duplicate removal\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen>(amplusplus::counter_coalesced_message_type_gen(1 << 12)));
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, simple cache duplicate removal (size 16)\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing>(amplusplus::counter_coalesced_message_type_gen(1 << 12), 4, amplusplus::no_routing(trans.rank(), trans.size())));
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, simple cache duplicate removal (size 512)\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing>(amplusplus::counter_coalesced_message_type_gen(1 << 12), 9, amplusplus::no_routing(trans.rank(), trans.size())));
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, simple cache duplicate removal (size 2048)\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing>(amplusplus::counter_coalesced_message_type_gen(1 << 12), 11, amplusplus::no_routing(trans.rank(), trans.size())));
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, per_thread cache duplicate removal (size 16)\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing>(amplusplus::counter_coalesced_message_type_gen(1 << 12), 4, amplusplus::no_routing(trans.rank(), trans.size())));
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, per_thread cache duplicate removal (size 512)\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing>(amplusplus::counter_coalesced_message_type_gen(1 << 12), 9, amplusplus::no_routing(trans.rank(), trans.size())));
  }

  if (rank == 0) fprintf(stderr, "One transport per iteration, per_thread cache duplicate removal (size 2048)\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, amplusplus::per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing>(amplusplus::counter_coalesced_message_type_gen(1 << 12), 11, amplusplus::no_routing(trans.rank(), trans.size())));
  }

#if 0
  if (rank == 0) fprintf(stderr, "One transport per iteration, unordered_set duplicate removal\n");
  for (int i = 0; i < 10; ++i) {
    amplusplus::transport trans = env.create_transport();
    do_test(my_graph, trans, my_size, my_start, my_end, chunk_size, graph_size, num_threads, reduction_unordered_set());
  }
#endif

  return 0;
}

int main(int argc, char** argv) {
#if IS_MPI_TRANSPORT
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  return do_one_thread(argc, argv, env);
#elif IS_SHM_TRANSPORT
  boost::scoped_ptr<amplusplus::shm_environment_common> common;

  int result = 0;
#pragma omp parallel firstprivate(result)
  {
#pragma omp single
    {
      common.reset(new amplusplus::shm_environment_common(omp_get_num_threads()));
      reduction_barrier.reset(new boost::barrier(omp_get_num_threads()));
    }
    amplusplus::environment env = amplusplus::shm_environment(*common, omp_get_thread_num());
    result = do_one_thread(argc, argv, env); // Arbitrary value is OK
  }
  return result;
#endif
}
