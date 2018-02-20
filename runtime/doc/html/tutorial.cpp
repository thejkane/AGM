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
#include <boost/property_map/property_map.hpp>

struct empty_deleter {
  typedef void result_type;
  void operator()() const {}
  template <typename T> void operator()(const T&) const {}
};

struct my_handler {
  // Raw handler
  void operator()(int source, const int* buf, int /*count*/) const {};
  // Coalesced handler
  void operator()(int source, int data) const {};
  // Object-based addressing
  void operator()(int data) const {};
  // Per-thread cache
  void operator()(std::pair<int,int>) const {};
};

struct owner_map_type {};

template<typename T>
amplusplus::transport::rank_type 
get(const owner_map_type, const T) { return 0; } // something more complicated

int main(int argc, char** argv) {
  (void)argc; (void)argv;

  using namespace amplusplus;
  using namespace boost;
  bool use_threads = true;
  owner_map_type owner_map;
  environment env = mpi_environment(argc, argv, use_threads);
  {
    // transport is reference-counted and thus copyable (using a shallow copy)
    transport trans = env.create_transport();
    // Get size and rank
    transport::rank_type r = trans.rank(), s = trans.size();
    // Non-coalesced message type
    message_type<int> mt = trans.create_message_type<int>();
    mt->set_max_count(100);
    mt.set_handler(my_handler());
    // Coalesced message type
    basic_coalesced_message_type<int, my_handler> mt2(basic_coalesced_message_type_gen(100), trans);
    mt2.set_handler(my_handler());
    // Object-based addressing
    typedef simple_generator<basic_coalesced_message_type_gen> simple_gen;
    simple_gen::
      call_result<int, my_handler, owner_map_type, no_reduction_t>::type
      mt3(simple_gen(basic_coalesced_message_type_gen(100)),
          trans, owner_map, no_reduction);
    // Routing
    typedef routing_generator<basic_coalesced_message_type_gen, ring_routing> routing_gen;
    routing_gen::
      call_result<int, my_handler, owner_map_type, no_reduction_t>::type
      mt4(routing_gen(basic_coalesced_message_type_gen(100), ring_routing(trans.rank(), trans.size())), trans, owner_map, no_reduction);
    // Caching, counter coalescing
    typedef per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> cache_gen;
    cache_gen::call_result<std::pair<int, int>, my_handler, owner_map_type, idempotent_combination_t<std::greater<int> > >::type mt5(cache_gen(counter_coalesced_message_type_gen(1024), 20, no_routing(trans.rank(), trans.size())), trans, owner_map, idempotent_combination(std::greater<int>()));
    {
      trans.set_nthreads(3); // Set number of threads that must enter and leave epoch
      transport::rank_type dest = 0;

      shared_array<int> buf(new int[3]);
      // Epoch control (explicit begin_epoch, end_epoch are deprecated)
      {
        scoped_epoch e(trans); // Begin epoch
        mt.message_being_built(0); // Mark that a message is coming; only needed without coalescing, but allowed elsewhere
        mt.send(buf.get(), 3, dest, empty_deleter() /* To control lifetime */);
        mt2.send(1, dest);
        mt3.send(2);
        mt4.send(3);
	mt5.send(std::make_pair(1,1));
      } // End epoch
      uintmax_t result = 0;
      {
        uintmax_t val;
        scoped_epoch_value e(trans, val, result); // Begin epoch
        mt.send(buf.get(), 3, dest, empty_deleter() /* To control lifetime */);
        val = 5;
      } // End epoch, accumulating into result

      // Old-style begin_epoch, end_epoch
      trans.begin_epoch();
      // ...
      trans.end_epoch();
    }
  } // Destroy transport
  return 0;
}
