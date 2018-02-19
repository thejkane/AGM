// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#ifndef FAST_RMAT_HPP
#define FAST_RMAT_HPP

#include <boost/graph/random/nway_bernoulli_distribution.hpp>
#include <boost/graph/random/splittable_ecuyer1988.hpp>
#include <boost/graph/random/btrd_binomial_distribution.hpp>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/graph/random/next_level_parameters_generated.ipp>
#include <stack>

// #define KEEP_ALL

//
// Available filters
//

namespace rmat_detail {

// keep all edges
template <typename T>
struct keep_all {
  bool operator()(T base_src, T base_tgt, T log2_nverts) const { return true; }
};

// keep only edges originating from a given source range
template <typename T>
struct source_range_filter {
  T start, end;
  source_range_filter(T start, T end): start(start), end(end) {}
  source_range_filter(): start(0), end(0) {}
  bool operator()(T base_src, T base_tgt, T log2_nverts) const {
    return base_src < end &&
      (base_src + (static_cast<T>(1) << log2_nverts)) > start &&
      base_src >= base_tgt &&
      (log2_nverts != 0 || base_src > base_tgt);  // Remove self-loops
  }
};

} // end namespace rmat_detail

// TODO: Replace this with boost::spirit::tuple?
template <typename T>
struct generated_edge {
  T src;
  T tgt;
  T multiplicity;
  generated_edge(T src, T tgt, T multiplicity = 1)
    : src(src), tgt(tgt), multiplicity(multiplicity) {}
  generated_edge(): src(0), tgt(0), multiplicity(0) {}
  bool operator==(const generated_edge<T>& e) const {
    return src == e.src && tgt == e.tgt && multiplicity == e.multiplicity;
  }
};

using boost::graph::random::nway_bernoulli_distribution;
using boost::graph::random::make_nway_bernoulli_distribution;

template <typename Gen>
static inline void alter_params(double& a, double& b, double& c, double& d, Gen* gen) {
  boost::variate_generator<Gen*, boost::uniform_01<> > gen01(gen, boost::uniform_01<>());
  a *= .9 + .2 * gen01();
  b *= .9 + .2 * gen01();
  c *= .9 + .2 * gen01();
  d *= .9 + .2 * gen01();
  double divisor = a + b + c + d;
  a /= divisor;
  b /= divisor;
  c /= divisor;
  d /= divisor;
}

template <typename Gen, typename T>
static inline
generated_edge<T> make_one_edge(T base_src, T base_tgt, int log2_nverts, Gen* gen, 
				double a, double b, double c, double d) {
  double params[4] = {a, b, c, d};
  while (log2_nverts != 0) {
    int quad = make_nway_bernoulli_distribution<4, double>(params, *gen)(*gen);
    --log2_nverts;
    T new_block_size = static_cast<T>(1) << log2_nverts;
    base_src += ((quad & 2) != 0 ? new_block_size : 0);
    base_tgt += ((quad & 1) != 0 ? new_block_size : 0);
    alter_params(a, b, c, d, gen);
  }
  return generated_edge<T>(base_src, base_tgt, 1);
}

template <typename Gen, typename T>
static inline
void make_quadrant_counts(T num_edges, T /* log2_nverts */, Gen* gen, double a, double b, 
			  double c, double d, T* quadrant_counts) {
  if (num_edges <= 20) {
    quadrant_counts[0] = 0;
    quadrant_counts[1] = 0;
    quadrant_counts[2] = 0;
    quadrant_counts[3] = 0;
    double params[4] = {a, b, c, d};

    nway_bernoulli_distribution<4, double, typename Gen::result_type> distrib 
      = make_nway_bernoulli_distribution<4, double>(params, *gen);

    for (T i = 0; i < num_edges; ++i) 
      ++quadrant_counts[distrib(*gen)];

  } else {

    boost::variate_generator<Gen*, boost::uniform_01<> > gen01(gen, boost::uniform_01<>());
    quadrant_counts[0] = btrd_binomial_distribution(num_edges, a, gen01);
    quadrant_counts[1] = btrd_binomial_distribution(num_edges - quadrant_counts[0], 
						    b / (1. - a), 
						    gen01);
    quadrant_counts[2] = btrd_binomial_distribution(num_edges - quadrant_counts[0] - quadrant_counts[1], 
						    c / (1. - a - b), 
						    gen01);
    quadrant_counts[3] = num_edges - quadrant_counts[0] - quadrant_counts[1] - quadrant_counts[2];
  }
}

namespace boost {

// Note: Generator must be a splittable generator such as the splittable_ecuyer1988
template <typename Transport, typename Generator, typename Graph, 
	  typename Filter = 
#ifdef KEEP_ALL
	  rmat_detail::keep_all<typename graph_traits<Graph>::vertices_size_type> >
#else
	  rmat_detail::source_range_filter<typename graph_traits<Graph>::vertices_size_type> >
#endif
class recursive_rmat_generator
  : public boost::iterator_facade<recursive_rmat_generator<Transport, Generator, Graph, Filter>, 
// 				  generated_edge<typename graph_traits<Graph>::vertices_size_type>, 
				  std::pair<typename graph_traits<Graph>::vertices_size_type,
					    typename graph_traits<Graph>::vertices_size_type>,
				  std::forward_iterator_tag, 
// 				  const generated_edge<typename graph_traits<Graph>::vertices_size_type>&>
				  const std::pair<typename graph_traits<Graph>::vertices_size_type,
						  typename graph_traits<Graph>::vertices_size_type>&>
{
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename graph_traits<Graph>::edges_size_type edges_size_type;

public:

  typedef std::pair<typename graph_traits<Graph>::vertices_size_type,
		    typename graph_traits<Graph>::vertices_size_type> value_type;

//   recursive_rmat_generator(const recursive_rmat_generator& /* x */) { assert(false); exit(-1); }
//   recursive_rmat_generator& operator=(const recursive_rmat_generator& /* x */) { assert(false); exit(-1); }

  recursive_rmat_generator(const Transport& transport, Generator& gen, vertices_size_type log2_nverts, 
			   vertices_size_type nedges, double a, double b, double c, double d) 
    : work(), filter(
#ifdef KEEP_ALL
		     rmat_detail::keep_all<vertices_size_type>()),
#else
                     rmat_detail::source_range_filter<vertices_size_type>
                       (transport.rank() * pow(2, log2_nverts) / transport.size(),
		       std::min((transport.rank() + 1) * pow(2, log2_nverts) / transport.size(), 
		       pow(2, log2_nverts)))), 
#endif
      param_data(4 * (log2_nverts + 1)), cached_edge(), at_end(false)
  {
    assert (a >= b);
    assert (a >= c);
    assert (b >= d);
    assert (c >= d);
    work_queue_entry s;
    s.gen = *gen.split_off_n(1).first;
    s.num_edges = nedges;
    s.log2_nverts = log2_nverts;
    s.base_src = s.base_tgt = 0;
    s.a = a;
    s.b = b;
    s.c = c;
    s.d = d;
    work.push(s);
    at_end = !step();
  }

  recursive_rmat_generator(const Transport& transport, Generator& gen, vertices_size_type log2_nverts, 
			   vertices_size_type nedges, double a, double b, double c, double d, 
			   const Filter& filter)
    : work(), filter(filter), param_data(4 * (log2_nverts + 1)),
      cached_edge(), at_end(false)
  {
    assert (a >= b);
    assert (a >= c);
    assert (b >= d);
    assert (c >= d);
    work_queue_entry s;
    s.gen = *gen.split_off_n(1).first;
    s.num_edges = nedges;
    s.log2_nverts = log2_nverts;
    s.base_src = s.base_tgt = 0;
    s.a = a;
    s.b = b;
    s.c = c;
    s.d = d;
    work.push(s);
    at_end = !step();
  }

  recursive_rmat_generator()
    : work(), filter(Filter()), param_data(), cached_edge(), at_end(true)
  { }

  // iterator_facade operations
  friend class boost::iterator_core_access;

private:

//   const generated_edge<vertices_size_type>& dereference() const {return cached_edge;}
  const value_type& dereference() const { return current_edge; }
  vertices_size_type get_multiplicity() const { return cached_edge.multiplicity; }

  bool equal(const recursive_rmat_generator& a) const {
    return at_end == a.at_end && work == a.work;
  }

  void increment() {
    if (at_end) return;
    at_end = !step();
  }

  // Internals
  struct work_queue_entry {
    Generator gen;
    edges_size_type num_edges;
    vertices_size_type log2_nverts;
    vertices_size_type base_src, base_tgt;
    double a, b, c, d;

    friend bool operator==(const work_queue_entry& a, const work_queue_entry& b) {
      return a.num_edges == b.num_edges &&
             a.log2_nverts == b.log2_nverts &&
             a.base_src == b.base_src && a.base_tgt == b.base_tgt;
    }
  };

  std::stack<work_queue_entry, std::vector<work_queue_entry> > work;
  Filter filter;
  std::vector<double> param_data;
  generated_edge<vertices_size_type> cached_edge;
  value_type current_edge;
  bool at_end;

  bool step() { // Returns true if there is a new edge
    while (true) { // Exited by return statements in loop
      if (work.empty()) return false;
      work_queue_entry e = work.top();
      work.pop();

      if (e.log2_nverts == 0) {
        assert (e.num_edges != 0);
        cached_edge = generated_edge<vertices_size_type>(e.base_src, e.base_tgt, e.num_edges);
	current_edge = std::make_pair(cached_edge.src, cached_edge.tgt);
        return true;
      } else if (e.num_edges == 1) {
        cached_edge = make_one_edge(e.base_src, e.base_tgt, e.log2_nverts, 
				    &e.gen, e.a, e.b, e.c, e.d);
	current_edge = std::make_pair(cached_edge.src, cached_edge.tgt);

        if (filter(cached_edge.src, cached_edge.tgt, 0)) {
          return true;
        } 
      } else {
        edges_size_type quadrant_counts[4];
        make_quadrant_counts(e.num_edges, e.log2_nverts, &e.gen, e.a, e.b, e.c, e.d, quadrant_counts);
        typename Generator::split_iterator subparts = e.gen.split_n(4).first;
        --e.log2_nverts;
        vertices_size_type old_base_src = e.base_src;
        vertices_size_type old_base_tgt = e.base_tgt;
        vertices_size_type new_block_size = static_cast<vertices_size_type>(1) << e.log2_nverts;
        double na = e.a, nb = e.b, nc = e.c, nd = e.d;

        for (int i = 3; i >= 0; --i) { // Order reversed to push in stack correctly
          e.gen = *subparts++;
          if (quadrant_counts[i] != 0) {
            e.num_edges = quadrant_counts[i];
            e.base_src = old_base_src + ((i & 2) ? new_block_size : 0);
            e.base_tgt = old_base_tgt + ((i & 1) ? new_block_size : 0);
            if (filter(e.base_src, e.base_tgt, e.log2_nverts)) {
              e.a = na;
              e.b = nb;
              e.c = nc;
              e.d = nd;
              alter_params(e.a, e.b, e.c, e.d, &e.gen);
              work.push(e);
            }
          }
        }
        // Go to top of loop
      }
    }
  }
};

} // end namespace boost

#endif // FAST_RMAT_HPP
