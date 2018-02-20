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

#ifndef BOOST_DISTRIBUTED_ERDOS_RENYI_GENERATOR_HPP
#define BOOST_DISTRIBUTED_ERDOS_RENYI_GENERATOR_HPP

#include "splittable_ecuyer1988.hpp"
#include <boost/iterator.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/uniform_01.hpp>

namespace boost {
  template <typename VertexDistribution, typename Generator = boost::random::splittable_ecuyer1988, typename VertexType = size_t>
  class distributed_erdos_renyi_iterator:
          public boost::iterator_facade<distributed_erdos_renyi_iterator<VertexDistribution, Generator, VertexType>,
                                        std::pair<VertexType, VertexType>,
                                        std::input_iterator_tag,
                                        std::pair<VertexType, VertexType> >
  {
    VertexDistribution owner;
    typedef size_t rank_type;
    rank_type my_rank;
    typename Generator::split_iterator_pair generators; // One for each source vertex
    Generator current_generator;
    VertexType current_source, current_target, nverts;
    bool is_undirected, allow_self_loops;
    boost::geometric_distribution<VertexType> distrib;

    void next_source_vertex() {
      // fprintf(stderr, "next_source_vertex()\n");
      ++current_source;
      current_target = (VertexType)(-1);
      ++generators.first;
      if (generators.first != generators.second) {
        current_generator = *generators.first;
      }
    }

    void step_to_owned_vertex() {
      // Find a source that I own
      while (current_source != nverts &&
             rank_type(owner(current_source)) != my_rank) {
        next_source_vertex();
      }
    }

    void raw_increment() { // Does not always produce a valid vertex
      if (current_source == nverts) return;
      boost::random::uniform_01_wrapper<Generator> wr(current_generator);
      VertexType delta = distrib(wr);
      // fprintf(stderr, "delta = %zu\n", delta);
      // First test is to deal with overflows, plus handle that the current
      // target might be -1 and so delta == nverts would be valid
      if (delta > nverts || current_target + delta >= nverts) {
        next_source_vertex();
        step_to_owned_vertex();
        // Return target as -1, so raw_increment() will be repeated
      } else {
        current_target += delta;
      }
      // fprintf(stderr, "New loc is (%zu, %zu)\n", current_source, current_target);
    }

    public:
    distributed_erdos_renyi_iterator(): owner(), my_rank(0), current_source(0), current_target((VertexType)(-1)), nverts(0) {}

    distributed_erdos_renyi_iterator(Generator& gen, VertexType nverts, double prob, bool is_undirected, bool allow_self_loops, VertexDistribution owner, rank_type my_rank): owner(owner), my_rank(my_rank), current_source(0), current_target((VertexType)(-1)), nverts(nverts), is_undirected(is_undirected), allow_self_loops(allow_self_loops), distrib(1. - prob) {
      generators = gen.split_off_n(nverts);
      raw_increment();
    }

    std::pair<VertexType, VertexType> dereference() const {
      return std::make_pair(current_source, current_target);
    }

    bool equal(const distributed_erdos_renyi_iterator& o) const {
      // Only really works on end iterators or those generated with the same parameters
      return (nverts - current_source) == (o.nverts - o.current_source) && current_target == o.current_target;
    }

    void increment() {
      do {
        raw_increment();
      } while (current_source != nverts &&
               (current_target == (VertexType)(-1) ||
                (!allow_self_loops && current_source == current_target) ||
                (is_undirected && current_target > current_source)));
    }
  };

}

#endif // BOOST_DISTRIBUTED_ERDOS_RENYI_GENERATOR_HPP
