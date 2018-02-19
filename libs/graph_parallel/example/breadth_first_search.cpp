// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

// Example usage of breadth_first_search algorithm

#include <am++/am++.hpp>

// Communicate via MPI
#include <am++/mpi_transport.hpp>

// Enable PBGL interfaces to BGL algorithms
#include <boost/graph/use_mpi.hpp>

// Breadth-first search algorithm
#include <boost/graph/breadth_first_search.hpp>

// Distributed adjacency list
#include <boost/graph/distributed/adjacency_list.hpp>

// METIS Input
#include <boost/graph/metis.hpp>

// Graphviz Output
#include <boost/graph/distributed/graphviz.hpp>

// Standard Library includes
#include <fstream>
#include <string>

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

typedef int weight_type;

// Edge property bundle
struct WeightedEdge {
  WeightedEdge(weight_type weight = 0) : weight(weight) { }
  
  weight_type weight;
};

// Vertex property bundle
struct VertexProperties {
  VertexProperties(int d = 0) : distance(d) { }

  int distance;
};

// Define MPI datatypes for bundled properties so AM++ can send messages containing them
namespace amplusplus {
  template<>
  struct make_mpi_datatype<WeightedEdge> : make_mpi_datatype_base {
    make_mpi_datatype<weight_type> dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1() {
      int blocklengths[1] = {1};
      MPI_Aint displacements[1];
      WeightedEdge test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.weight, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Datatype types[1] = {dt1.get()};
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };

  template<>
  struct make_mpi_datatype<VertexProperties> : make_mpi_datatype_base {
    make_mpi_datatype<int> dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1() {
      int blocklengths[1] = {1};
      MPI_Aint displacements[1];
      VertexProperties test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.distance, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Datatype types[1] = {dt1.get()};
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}


using namespace boost;

/* An undirected graph with distance values stored on the vertices. */
typedef adjacency_list<vecS, distributedS<vecS>, directedS,
                       VertexProperties, WeightedEdge>
  Graph;

int main(int argc, char* argv[])
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  amplusplus::transport trans = env.create_transport();

  // Parse command-line options
  const char* filename = "weighted_graph.gr";
  if (argc > 1) filename = argv[1];

  // Open the METIS input file
  std::ifstream in(filename);
  graph::metis_reader reader(in);

  // Load the graph using the default distribution
  Graph g(reader.begin(), reader.end(), reader.num_vertices(), trans);

  // Get vertex 0 in the graph
  graph_traits<Graph>::vertex_descriptor start = vertex(0, g);

  // Compute BFS levels from vertex 0
  typedef property_map<Graph, int VertexProperties::*>::type DistanceMap;
  DistanceMap distance = get(&VertexProperties::distance, g);

  {
    amplusplus::scoped_epoch epoch(trans);
    put(distance, start, 0);
  }
  using boost::graph::distributed::breadth_first_search;
  breadth_first_search<Graph, bfs_visitor<distance_recorder<DistanceMap, on_tree_edge> > >
    bfs(g, make_bfs_visitor(record_distances(distance, on_tree_edge())));

  // Run bfs from start vertex
  bfs.run(start);

  // Output a Graphviz DOT file
  std::string outfile;

  if (argc > 2)
    outfile = argv[2];
  else {
    outfile = filename;
    size_t i = outfile.rfind('.');
    if (i > 0 && i != std::string::npos)
      outfile.erase(outfile.begin() + i, outfile.end());
    outfile += "-bfs.dot";
  }

  if (trans.rank() == 0) {
    std::cout << "Writing GraphViz output to " << outfile << "... ";
    std::cout.flush();
  }
  write_graphviz(outfile, g,
                 make_label_writer(distance));
  if (trans.rank() == 0)
    std::cout << "Done." << std::endl;

  return 0;
}
