// Copyright (C) 2009-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Andrew Lumsdaine

#include <am++/am++.hpp> 
#include <am++/mpi_transport.hpp>

// A test of the distributed compressed sparse row graph type
#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/test/minimal.hpp>

/****************************************************************************
 * Edge weight generator iterator                                           *
 ****************************************************************************/
template<typename F, typename RandomGenerator>
class generator_iterator
{
public:
  typedef std::input_iterator_tag iterator_category;
  typedef typename F::result_type value_type;
  typedef const value_type&       reference;
  typedef const value_type*       pointer;
  typedef void                    difference_type;

  explicit 
  generator_iterator(RandomGenerator& gen, const F& f = F()) 
    : f(f), gen(&gen) 
  { 
    value = this->f(gen); 
  }

  reference operator*() const  { return value; }
  pointer   operator->() const { return &value; }

  generator_iterator& operator++()
  {
    value = f(*gen);
    return *this;
  }

  generator_iterator operator++(int)
  {
    generator_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const generator_iterator& other) const
  { return f == other.f; }

  bool operator!=(const generator_iterator& other) const
  { return !(*this == other); }

private:
  F f;
  RandomGenerator* gen;
  value_type value;
};

template<typename F, typename RandomGenerator>
inline generator_iterator<F, RandomGenerator> 
make_generator_iterator( RandomGenerator& gen, const F& f)
{ return generator_iterator<F, RandomGenerator>(gen, f); }

using namespace boost;

typedef int weight_type;

struct EdgeProperty {
  EdgeProperty(weight_type weight = 0) : weight(weight) { }
  
  weight_type weight;
};

struct VertexProperty {
  VertexProperty(int d = 0)
    : distance(d) { }

  int distance;
};

int test_main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(1);

  int n = 40;
  double prob = 0.1;
  int c = 10;

  minstd_rand gen;

  // No properties
  {
    typedef compressed_sparse_row_graph<directedS, no_property, no_property, 
                                        no_property, distributedS<> >
      Digraph;

    parallel::block<graph_traits<Digraph>::vertices_size_type> dist(trans, n);

    typedef sorted_erdos_renyi_iterator<minstd_rand, Digraph> SortedERIter;
    typedef erdos_renyi_iterator<minstd_rand, Digraph> ERIter;

    std::vector<ERIter::value_type> edges;
    std::vector<graph_traits<Digraph>::vertex_descriptor> sources, targets;
    ERIter start(gen, n, prob);
    for ( ; start != ERIter() ; ++start) {
      edges.push_back(*start);
      
      if (dist(start->first) == trans.rank()) { 
        sources.push_back(start->first);
        targets.push_back(start->second);
      }
    }

    Digraph g(edges_are_unsorted, ERIter(gen, n, prob), ERIter(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs(g);
    bfs.run(vertex(0, g));

    Digraph g2(edges_are_sorted, SortedERIter(gen, n, prob), SortedERIter(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs2(g2);
    bfs2.run(vertex(0, g2));

    Digraph g3(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs3(g3);
    bfs3.run(vertex(0, g3));

    Digraph g4(distributed_construct_inplace_from_sources_and_targets, sources, targets, n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs4(g4);
    bfs4.run(vertex(0, g4));

    Digraph g5(n, trans, parallel::oned_block_cyclic<graph_traits<Digraph>::vertices_size_type>(trans, 1));
    add_edges(edges.begin(), edges.end(), g5);
    boost::graph::distributed::breadth_first_search<Digraph> bfs5(g5);
    bfs5.run(vertex(0, g5));

    // TODO (NGE): Test add_edges_sorted()?
  }

  // Edge properties
  {
    typedef compressed_sparse_row_graph<directedS, no_property, EdgeProperty, no_property, 
                                        distributedS<> >
      Digraph;

    parallel::block<graph_traits<Digraph>::vertices_size_type> dist(trans, n);

    typedef sorted_erdos_renyi_iterator<minstd_rand, Digraph> SortedERIter;
    typedef erdos_renyi_iterator<minstd_rand, Digraph> ERIter;

    std::vector<ERIter::value_type> edges;
    std::vector<graph_traits<Digraph>::vertex_descriptor> sources, targets;
    std::vector<EdgeProperty> edge_props;
    uniform_int<int> edge_weight(0, c);

    ERIter start(gen, n, prob);
    for ( ; start != ERIter() ; ++start) {
      edges.push_back(*start);
      
      if (dist(start->first) == trans.rank()) { 
        sources.push_back(start->first);
        targets.push_back(start->second);
        edge_props.push_back(EdgeProperty(edge_weight(gen)));
      }
    }

    Digraph g(edges_are_unsorted, ERIter(gen, n, prob), ERIter(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs(g);
    bfs.run(vertex(0, g));

    Digraph g2(edges_are_sorted, SortedERIter(gen, n, prob), SortedERIter(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs2(g2);
    bfs2.run(vertex(0, g2));

    Digraph g3(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), 
	       make_generator_iterator(gen, uniform_int<int>(0, c)), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs3(g3);
    bfs3.run(vertex(0, g3));

    Digraph g4(distributed_construct_inplace_from_sources_and_targets, sources, 
	       targets, edge_props, n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs4(g4);
    bfs4.run(vertex(0, g4));

    Digraph g5(n, trans, parallel::oned_block_cyclic<graph_traits<Digraph>::vertices_size_type>(trans, 1));
    add_edges(edges.begin(), edges.end(), g5);
    boost::graph::distributed::breadth_first_search<Digraph> bfs5(g5);
    bfs5.run(vertex(0, g5));

    // TODO (NGE): Test add_edges_sorted()?
  }

  // Edge and vertex properties
  {
    typedef compressed_sparse_row_graph<directedS, VertexProperty, EdgeProperty, 
                                        no_property, 
                                        distributedS<> >
      Digraph;

    parallel::block<graph_traits<Digraph>::vertices_size_type> dist(trans, n);

    typedef sorted_erdos_renyi_iterator<minstd_rand, Digraph> SortedERIter;
    typedef erdos_renyi_iterator<minstd_rand, Digraph> ERIter;

    std::vector<ERIter::value_type> edges;
    std::vector<graph_traits<Digraph>::vertex_descriptor> sources, targets;
    std::vector<EdgeProperty> edge_props;
    uniform_int<int> edge_weight(0, c);

    ERIter start(gen, n, prob);
    for ( ; start != ERIter() ; ++start) {
      edges.push_back(*start);
      
      if (dist(start->first) == trans.rank()) { 
        sources.push_back(start->first);
        targets.push_back(start->second);
        edge_props.push_back(EdgeProperty(edge_weight(gen)));
      }
    }

    Digraph g(edges_are_unsorted, ERIter(gen, n, prob), ERIter(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs(g);
    bfs.run(vertex(0, g));

    Digraph g2(edges_are_sorted, SortedERIter(gen, n, prob), SortedERIter(), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs2(g2);
    bfs2.run(vertex(0, g2));

    Digraph g3(edges_are_unsorted_multi_pass, edges.begin(), edges.end(), 
	       make_generator_iterator(gen, uniform_int<int>(0, c)), n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs3(g3);
    bfs3.run(vertex(0, g3));


    Digraph g4(distributed_construct_inplace_from_sources_and_targets, sources, 
               targets, edge_props, n, trans, dist);
    boost::graph::distributed::breadth_first_search<Digraph> bfs4(g4);
    bfs4.run(vertex(0, g4));

    // NGE: Both add_edges() and bfs5.run() below fail
#if 0 
    Digraph g5(n, trans, parallel::oned_block_cyclic(trans, 1));
    add_edges(edges.begin(), edges.end(), edge_props.begin(), edge_props.end(), g5);
    boost::graph::distributed::breadth_first_search<Digraph> bfs5(g5);
    bfs5.run(vertex(0, g5));
#endif

    // TODO (NGE): Test add_edges_sorted()?
  }

  return 0;
}
