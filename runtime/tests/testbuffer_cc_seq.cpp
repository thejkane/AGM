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
//  Authors: Thejaka Kanewala
//           Jeremiah Willcock
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

typedef uint64_t Arg_t;
#define MAX_T_SIZE 1000

class variable_buffer_test {
public:

  void receive(Arg_t a, int tid) {
    ++receive_count;
  }

  struct buffer_handler {
  private:
    variable_buffer_test* bftest;
  public:
    buffer_handler() : bftest(NULL) {}
    buffer_handler(variable_buffer_test* bft) : bftest(bft) {}
    void operator()(amplusplus::transport::rank_type& src, const Arg_t& arg) {
      assert(src == 0);
      int tid = amplusplus::detail::get_thread_id();
      bftest->receive(arg, tid);
    }
  };

  typedef amplusplus::counter_coalesced_message_type_gen::inner<Arg_t, buffer_handler>::type msg_type;

  amplusplus::transport trans;
  msg_type mt;
  Arg_t receive_count;

  
  variable_buffer_test(amplusplus::transport& t):trans(t),
						 mt(amplusplus::counter_coalesced_message_type_gen(1 << 12), t), receive_count(0){
  }


  void initialize() {
    mt.set_handler(buffer_handler(this));
  }

  void inline operator()(int tid) {
    run(tid);
  }

  void run(int tid) {
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      amplusplus::rank_type dest = 1;
      { 
	amplusplus::scoped_epoch epoch(trans); 
	if (trans.rank() == 0) {
	  for (int i=0; i < MAX_T_SIZE; ++i) {
	    mt.send((Arg_t)i, dest);
	  }
	}
      }

      std::cout << "Receive count = " << receive_count << std::endl;
    
      //int blocksz = 100;

      //Vertex* buff = (Vertex*)malloc(sizeof(Vertex) * (blocksz+1));

      //buff[0] = 107;
      //for (int i=0; i < 100; ++i) {
      //  buff[i+1] = i;
      //}

      //amplusplus::rank_type r = tm.get_transport().rank();
      //amplusplus::rank_type dest;

      //dest = r % 2;
      //tm.send(buff, (blocksz+1), dest, boost::bind(empty_deleter(), msg));
    }
  }



};


int main(int argc, char** argv) {
  (void)argc; (void)argv;
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  amplusplus::transport trans = env.create_transport();
  variable_buffer_test vbt(trans);
  vbt.initialize();
  vbt(0);
  std::cout << "End of test ..." << std::endl;
  return 0;
}



