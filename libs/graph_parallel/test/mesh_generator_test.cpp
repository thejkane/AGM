// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/graph/mesh_graph_generator.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <boost/test/minimal.hpp>

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

using namespace boost;

/****************************************************************************
 * Timing
 ****************************************************************************/
typedef double time_type;

inline time_type get_time()
{
        return MPI_Wtime();
}

std::string print_time(time_type t)
{
  std::ostringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
  return out.str();
}

/****************************************************************************
 * Edge weight generator iterator                                           *
 ****************************************************************************/
template<typename F>
class generator_iterator
{
public:
  typedef std::input_iterator_tag iterator_category;
  typedef typename F::result_type value_type;
  typedef const value_type&       reference;
  typedef const value_type*       pointer;
  typedef void                    difference_type;

  explicit generator_iterator(const F& f = F()) : f(f) { value = this->f(); }

  reference operator*() const  { return value; }
  pointer   operator->() const { return &value; }

  generator_iterator& operator++()
  {
    value = f();
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
  value_type value;
};

template<typename F>
inline generator_iterator<F> make_generator_iterator(const F& f)
{ return generator_iterator<F>(f); }

int test_main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(1);

  if (argc < 5) {
    std::cerr << "Usage: mesh_generator_test <x> <y> <toroidal> <emit dot file>\n";
    exit(-1);
  }

  int x(atoi(argv[1])), y(atoi(argv[2]));
  bool toroidal(argv[3] == std::string("true"));
  bool emit_dot_file(argv[4] == std::string("true"));

  typedef adjacency_list<vecS, 
                         distributedS<vecS>,
                         directedS> Graph;

  Graph g(mesh_iterator<Graph>(x, y, toroidal),
          mesh_iterator<Graph>(),
          x*y, trans);

  BGL_FORALL_VERTICES(v, g, Graph)
    if (toroidal) 
      assert(out_degree(v, g) == 4);
    else
      assert(out_degree(v, g) >= 2 && out_degree(v, g) <= 4);  

  if ( emit_dot_file )
    write_graphviz("mesh.dot", g);

  return 0;
}
