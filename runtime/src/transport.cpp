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

#include <am++/transport.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/move/move.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/property_maps/null_property_map.hpp>
#include <algorithm>
#include <set>
#include <map>

namespace amplusplus {

environment::environment(boost::shared_ptr<environment_base> env)
  : env(env)
{
}

environment::environment(const environment& e)
  : env(e.env)
{
}

transport environment::create_transport() {return env->create_transport(*this);}

template <typename E>
struct first_is_t {
  E val;
  explicit first_is_t(const E& val): val(val) {}
  template <typename P> bool operator()(const P& p) const {return p.first == val;}
};

template <typename E> first_is_t<E> first_is(const E& e) {return first_is_t<E>(e);}

void transport_base::add_flush_object(const boost::function<bool ()>& f) {
  boost::lock_guard<amplusplus::detail::mutex> l(flush_lock);
  flushes_needed.push_back(f);
}

scheduler::task_result transport_base::flush() {
  bool did_anything = false;
  std::vector<boost::function<bool ()> > tmp;
  while(true) {
    boost::function<bool ()> f;
    {
      boost::lock_guard<amplusplus::detail::mutex> l(flush_lock);
      if (flushes_needed.empty()) break;
      f.swap(flushes_needed.back());
      flushes_needed.pop_back();
    }
    assert (f);
    did_anything = true;
    if(f()) {
      tmp.push_back(f);
    }
  }
  if(!tmp.empty()) {
    boost::lock_guard<amplusplus::detail::mutex> l(flush_lock);
    flushes_needed.insert(flushes_needed.end(), tmp.begin(), tmp.end());
  }
  return did_anything ? scheduler::tr_busy : scheduler::tr_idle;
}

bool transport::operator==(const transport& o) const {return trans_base == o.trans_base;}
bool transport::operator!=(const transport& o) const {return trans_base != o.trans_base;}

}
