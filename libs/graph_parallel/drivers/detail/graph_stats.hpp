// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

#ifndef PBGL2_GRAPH_STATS
#define PBGL2_GRAPH_STATS

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <boost/graph/distributed/owner_defs.hpp>

#include <iostream>
#include <fstream>
#include <mpi.h>

template <typename Graph>
class graph_stats {

private:
  int rank;

public:
  graph_stats() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }

  void print_degree_distribution(const Graph& g) {
    std::string fname = "degree_dist_";
    fname += std::to_string(rank);
    fname += ".txt";
    
    std::ofstream degrfile;
    degrfile.open(fname);
    degrfile << "Vertices\t Degree\n";
    
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      degrfile << v << "\t" << boost::out_degree(v, g) << "\n";
    }

    degrfile.close();    
  }

};



#endif
