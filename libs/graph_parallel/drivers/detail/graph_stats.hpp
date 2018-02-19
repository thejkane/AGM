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
