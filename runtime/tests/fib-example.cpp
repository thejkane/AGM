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
//  Authors: Marcin Zalewski
//           Jeremiah Willcock
//           Andrew Lumsdaine


// This is a contrived example showing how fibonacci can be computed in AM++. Run it with something like:
//  openmpirun -np 2 -bynode -mca btl tcp,self ./tests/fib-example 2 28


#include <config.h>

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>
#include <am++/basic_coalesced_message_type.hpp>
#include <am++/detail/append_buffer.hpp>
#include <utility>
#include <vector>
#include <iostream>
#include <boost/thread/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <random>

struct fib_helper {
  fib_helper(bool is_root = false) : answer(0), is_root(is_root) {}
  fib_helper(size_t idx, amplusplus::transport::rank_type src) : answer(0), idx(idx), src(src), is_root(false) {
  }
  fib_helper(const fib_helper& f) : idx(f.idx), src(f.src), is_root(f.is_root)  {
    answer = f.answer.load();
  }
  fib_helper& operator=(const fib_helper& f) {
    idx = f.idx;
    src = f.src;
    is_root = f.is_root;
    answer = f.answer.load();
    return *this;
  }

  amplusplus::detail::atomic<unsigned int> answer;
  size_t idx;
  amplusplus::transport::rank_type src;
  bool is_root;
};

typedef std::pair<unsigned int, size_t> fib_data;

struct empty_deleter {
  typedef void result_type;
  void operator()() const {}
  template <typename T> void operator()(const T&) const {}
};

class fib {
public:
  fib(amplusplus::transport& trans, unsigned int nthreads);

  void operator()(unsigned int n, unsigned int tid) {
    assert(n > 2);
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      amplusplus::scoped_epoch epoch(trans);
      if(n % trans.size() == trans.rank() && tid == 0) {
	std::cout << "Beginning fib with " << n << std::endl;
	const size_t idx = responses.push_back(fib_helper(true));
	fib_message.send(std::make_pair(n-1, idx), (n-1)%trans.size());
	fib_message.send(std::make_pair(n-2, idx), (n-2)%trans.size());
      }
    }
  }

  void operator()(unsigned int tid) {
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      amplusplus::scoped_epoch epoch(trans);
    }
  }  

  void divide(const fib_data& data, const amplusplus::transport::rank_type src) {
    const size_t idx = responses.push_back(fib_helper(data.second, src));
    fib_message.send(std::make_pair(data.first-1, idx), distribution(generator));
    fib_message.send(std::make_pair(data.first-2, idx), distribution(generator));
  }

  void merge(const fib_data& data) {
    const unsigned int idx = data.second;
    while(responses[idx].answer.load() == 0) {
      unsigned int zero = 0;
      if(responses[idx].answer.compare_exchange_weak(zero, data.first))
	return;
    } 
    const unsigned int newval = responses[idx].answer.load() + data.first;
    if(responses[idx].is_root) {
      std::cout << "Answer: " << newval << std::endl;
      return;
    }
    response_message.send(std::make_pair(newval, responses[idx].idx), responses[idx].src);
  }

private:
  struct fib_handler;
  struct response_handler;

  amplusplus::transport& trans;
  amplusplus::basic_coalesced_message_type<fib_data, fib_handler> fib_message;
  amplusplus::basic_coalesced_message_type<fib_data, response_handler> response_message;
  amplusplus::detail::append_buffer<fib_helper>  responses;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution;
};

struct fib::fib_handler {
  fib_handler(fib& self) : self(&self) {}

  void operator()(const amplusplus::transport::rank_type src, const fib_data& data) const {
    if(data.first == 1 || data.first == 2) {
      self->response_message.send(std::make_pair(1, data.second), src);
    } else {
      self->divide(data, src);
    }
  }

  fib* self;
};

struct fib::response_handler {
  response_handler(fib& self) : self(&self) {}

  void operator()(const amplusplus::transport::rank_type src, const fib_data& data) const {
    self->merge(data);
  }

  fib* self;
};

fib::fib(amplusplus::transport& trans, unsigned int nthreads) : trans(trans), fib_message(amplusplus::basic_coalesced_message_type_gen(5), trans), response_message(amplusplus::basic_coalesced_message_type_gen(5), trans), responses(), generator(), distribution(0, trans.size() - 1) {
    fib_message.set_handler(fib_handler(*this));
    response_message.set_handler(response_handler(*this));
  }

int main(int argc, char* argv[]) {
  const unsigned int nthreads = boost::lexical_cast<unsigned int>(argv[1]);
  const unsigned int input = boost::lexical_cast<unsigned int>(argv[2]);

  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(nthreads);

  amplusplus::register_mpi_datatype<fib_data>();

  fib f(trans, nthreads);

  boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
  for (int i = 0; i < nthreads - 1; ++i) {
    boost::thread thr(boost::ref(f), i + 1);
    threads[i].swap(thr);
  }

  f(input, 0);

  for (int i = 0; i < nthreads - 1; ++i)
    threads[i].join();
  
  return 0;
}
