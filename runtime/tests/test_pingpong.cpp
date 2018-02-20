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
#include <boost/config.hpp>
#define TRANSPORT_HEADER <am++/BOOST_JOIN(TRANSPORT, _transport).hpp>
#include TRANSPORT_HEADER
#include <am++/basic_coalesced_message_type.hpp>
#include <am++/counter_coalesced_message_type.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assert.hpp>
#include <stdio.h>
#include <string>
#include <sstream>
#include <numeric>
#include <unistd.h>

class timer {
  double start;
  std::string label;
  int reps;
  size_t size;
  bool print;

  public:
  timer(const std::string& label, int reps, size_t size, bool print): start(amplusplus::get_time()), label(label), reps(reps), size(size), print(print) {
    std::cerr << "Starting time " << label << std::endl << std::flush;
  }
  ~timer() {
    double end = amplusplus::get_time();
    if (!print) return;
    double time = (end - start) / reps;
    std::ostringstream os;
    os << label << ": " << time * 1.e6 << " us";
    if (size != 0) {
      os << " = " << (size / time / 1048576.) << " MiB/s";
    }
    os << std::endl;
    std::cerr << os.str() << std::flush;
  }
};

struct empty_deleter {
  typedef void result_type;
  void operator()() const {}
  template <typename T> void operator()(const T&) const {}
};

struct pingpong_handler_empty_reply {
  int reps;
  typedef amplusplus::message_type<int> tm_type;
  tm_type& tm;
  int mutable counter;
  pingpong_handler_empty_reply(int reps, tm_type& tm): reps(reps), tm(tm), counter(0) { std::cerr << "Constructor" << std::endl; }
  void operator()(int source, const int* /*buf*/, int count) const {
    // std::cerr << "pingpong_handler_empty_reply(" << reps << " of " << counter << " count=" << count << ")\n" << std::flush;
    if (source != 0) { // rank 0
      if(++counter != reps) {
      // std::cerr << "rank 0, counter = " << counter << std::endl;
	tm.message_being_built(1);
	tm.send(0, 0, 1, empty_deleter());
      }
    } else { // rank 1
      // std::cerr << "rank 1 counter = " << counter << std::endl;
      tm.message_being_built(0);
      tm.send(0, 0, 0, empty_deleter());
      ++counter;
    }
  }
};

struct pingpong_handler_non_coalesced {
  int reps;
  typedef amplusplus::message_type<int> tm_type;
  tm_type& tm;
  mutable int msg;
  pingpong_handler_non_coalesced(int reps, tm_type& tm): reps(reps), tm(tm) {}
  void operator()(int source, const int* buf, int /*count*/) const {
    int data = *buf;
    if (source != 0) { // rank 0
      if (data + 1 != reps) {
        // boost::shared_ptr<int> msg = boost::make_shared<int>(data + 1);
        msg = data + 1;
        tm.message_being_built(1);
        tm.send(&msg, 1, 1, empty_deleter());
      } 
    } else { // rank 1
      // boost::shared_ptr<int> msg = boost::make_shared<int>(data);
      msg = data;
      tm.message_being_built(0);
      tm.send(&msg, 1, 0, empty_deleter());
    }
  }
};

template <typename Pool>
struct pool_deleter {
  Pool& pool;
  pool_deleter(Pool& pool): pool(pool) {}
  void operator()(void* p) {
    // fprintf(stderr, "%d: pool free %p\n", gasnet_mynode(), p);
    pool.free(p);
  }
};

struct dummy_handler {
  void operator()(int /*source*/, const int* /*buf*/, int /*count*/) {}
};

struct pingpong_handler_non_coalesced_sized {
  int reps;
  typedef amplusplus::message_type<int> tm_type;
  tm_type& tm;
  size_t size;
  boost::shared_ptr<int> msg;

  pingpong_handler_non_coalesced_sized(int reps, tm_type& tm, size_t size): reps(reps), tm(tm), size(size)
    { msg = this->allocate(tm.get_transport(), size); }

  static boost::shared_ptr<int> allocate(const amplusplus::transport& transport, size_t size) {
    // boost::shared_ptr<int> p = boost::static_pointer_cast<int>(transport.alloc_memory(sizeof(int) * size));
    boost::shared_ptr<int> p(new int[size], boost::checked_array_deleter<int>());
    assert (p.get());
    return p;
  }

  void operator()(int source, const int* buf, int /*count*/) const {
    int data = *buf;
    // fprintf(stderr, "source = %d, data = %d\n", source, data);
    if (source != 0) { // rank 0
      if (data + 1 != reps) {
        *msg = data + 1;
        tm.message_being_built(1);
        tm.send(msg.get(), size, 1, boost::bind(empty_deleter(), msg));
      } 
    } else { // rank 1
      *msg = data;
      tm.message_being_built(0);
      tm.send(msg.get(), size, 0, boost::bind(empty_deleter(), msg));
    }
  }
};

struct pingpong_handler_coalesced_basic {
  int reps;
  typedef amplusplus::basic_coalesced_message_type<int, pingpong_handler_coalesced_basic> tm_type;
  tm_type* tm;
  pingpong_handler_coalesced_basic(): reps(0), tm(0) {}
  pingpong_handler_coalesced_basic(int reps, tm_type& tm): reps(reps), tm(&tm) {}
  void operator()(int source, int data) {
    if (source != 0) { // rank 0
      if (data + 1 != reps) {
        (*tm).message_being_built(1);
        (*tm).send(data + 1, 1);
      }
    } else { // rank 1
      (*tm).message_being_built(0);
      (*tm).send(data, 0);
    }
  }
};

void do_one_thread(amplusplus::environment& env) {
#if 0
  MPI_Init(0, 0);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (size != 2) {std::cerr << "Must have exactly two ranks" << std::endl; return 1;}
#endif

  // const int reps = 100000;
  const int reps = 100;

#if 0
  if (1) {
    timer t("MPI 1", reps * 2, sizeof(int), (rank == 0));
    int x = 1;
    for (int i = 0; i < reps; ++i) {
      if (rank == 0) {
        MPI_Send(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }
    }
  }

  if (1) {
    timer t("MPI 2", reps * 2, sizeof(int), (rank == 0));
    int x = 1;
    for (int i = 0; i < reps; ++i) {
      if (rank == 0) {
        MPI_Send(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }
    }
  }
#endif

  if (1) {
    amplusplus::transport trans = env.create_transport();
    BOOST_ASSERT (trans.size() >= 2);
    amplusplus::message_type<int> tm = trans.create_message_type<int>();
    tm.set_max_count(0);
    tm.set_handler(pingpong_handler_empty_reply(reps, tm));
    timer t("AM++ empty messages", reps * 2, 0, (trans.rank() == 0));
    {
      amplusplus::scoped_epoch epoch(trans);
      std::cerr << "Started empty on rank " << trans.rank() << std::endl;
      if (trans.rank() == 0) {
        tm.message_being_built(1);
        tm.send(0, 0, 1, empty_deleter());
	std::cerr << "Sent" << std::endl;
      }
    }
  }

  std::cerr << "Epoch done" << std::endl;

  if (1) {
    amplusplus::transport trans = env.create_transport();
    BOOST_ASSERT (trans.size() >= 2);
    amplusplus::message_type<int> tm = trans.create_message_type<int>();
    tm.set_max_count(1);
    tm.set_handler(pingpong_handler_non_coalesced(reps, tm));
    amplusplus::scoped_epoch epoch(trans);
    int msg = 0;
    {
      timer t("AM++ no coalescing", reps * 2, sizeof(int), (trans.rank() == 0));
      if (trans.rank() == 0) {
        // boost::shared_ptr<int> msg = boost::make_shared<int>(0);
        tm.message_being_built(1);
        tm.send(&msg, 1, 1, empty_deleter());
      }
    }
  }

  if (1) {
    // fprintf(stderr, "max_medium = %zu, max_long = %zu\n", (size_t)gasnet_AMMaxMedium(), (size_t)gasnet_AMMaxLongRequest());
    static const size_t sizes[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152}; // Measured in sizeof(int)
    for (size_t size_idx = 0; size_idx < sizeof(sizes) / sizeof(size_t); ++size_idx) {
      const size_t size = sizes[size_idx];
      // if (size * sizeof(int) > gasnet_AMMaxLongRequest()) continue;
      size_t actual_reps = reps;
      if (size >= 4096) actual_reps /= (size / 4096);
      if (actual_reps <= 1) actual_reps = 2;
      amplusplus::transport trans = env.create_transport();
      BOOST_ASSERT (trans.size() >= 2);
      amplusplus::message_type<int> tm = trans.create_message_type<int>();
      tm.set_max_count(size);
      tm.set_handler(pingpong_handler_non_coalesced_sized(actual_reps, tm, size));
      amplusplus::scoped_epoch epoch(trans);
      boost::shared_ptr<int> msg(pingpong_handler_non_coalesced_sized::allocate(tm.get_transport(), size));
      {
        timer t("AM++ no coalescing (size = " + boost::lexical_cast<std::string>(size * sizeof(int)) + ")", actual_reps * 2, size * sizeof(int), (trans.rank() == 0));
        if (trans.rank() == 0) {
          // boost::shared_ptr<int> msg = boost::make_shared<int>(0);
          *msg = 0;
          tm.message_being_built(1);
          tm.send(msg.get(), size, 1, boost::bind(empty_deleter(), msg));
        }
      }
    }
  }

  if (1) {
    amplusplus::transport trans = env.create_transport();
    BOOST_ASSERT (trans.size() >= 2);
    amplusplus::basic_coalesced_message_type<int, pingpong_handler_coalesced_basic> tm(amplusplus::basic_coalesced_message_type_gen(1), trans);
    tm.set_handler(pingpong_handler_coalesced_basic(reps, tm));
    amplusplus::scoped_epoch epoch(trans);
    {
      timer t("AM++ basic_coalesced_message_type", reps * 2, sizeof(int), (trans.rank() == 0));
      if (trans.rank() == 0) {
        tm.message_being_built(1);
        tm.send(0, 1);
      }
    }
  }

  if (1) {
    amplusplus::transport trans = env.create_transport();
    BOOST_ASSERT (trans.size() >= 2);
    amplusplus::basic_coalesced_message_type<int, pingpong_handler_coalesced_basic> tm(amplusplus::basic_coalesced_message_type_gen(1), trans);
    tm.set_handler(pingpong_handler_coalesced_basic(reps, tm));
    {
      timer t("AM++ basic_coalesced_message_type (using TD)", reps * 2, sizeof(int), (trans.rank() == 0));
      {
        amplusplus::scoped_epoch epoch(trans);
        if (trans.rank() == 0) {
          tm.message_being_built(1);
          tm.send(0, 1);
        }
      }
    }
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
