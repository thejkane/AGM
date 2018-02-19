// My first test
//


#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <iostream>
#include <string>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/parallel/append_buffer.hpp>

#include <am++/counter_coalesced_message_type.hpp>

//=======
#define PBGL_ACCOUNTING

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/adjacency_list.hpp> // must precede delta_stepping...
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <boost/parallel/append_buffer.hpp>

#include <am++/counter_coalesced_message_type.hpp>


//========


using namespace boost;

typedef int32_t weight_type; 

struct WeightedEdge {
  WeightedEdge(weight_type weight = 0) : weight(weight) { }
  
  weight_type weight;
};

struct VertexProperties {
  VertexProperties(int d = 0) : distance(d) { }

  weight_type distance;
};

//================== MPT Structures ===================
// TODO what are we exactly doing here ?
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
    make_mpi_datatype<weight_type> dt1;
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


//==============================================

//typedef adjacency_list<vecS, vecS, directedS, no_property, MyDistance, no_property, listS> Graph; 
typedef adjacency_list<vecS, distributedS<>, directedS, VertexProperties, WeightedEdge> Graph; 

void run_shortest_path_algorithm(amplusplus::transport trans, Graph g)
{
    // 0th vertex is the source
    Graph::vertex_descriptor source = vertex(0, g);

    // If g is connected run delta-stepping and verify
    boost::graph::distributed::delta_stepping_shortest_paths<Graph, 
      typename boost::property_map<Graph, weight_type VertexProperties::*>::type, 
      typename boost::property_map<Graph, weight_type WeightedEdge::*>::type> 
          D(g, get(&VertexProperties::distance, g), get(&WeightedEdge::weight, g));

    D.set_source(source);

    int num_threads = 5;
    // Many threads now
    trans.set_nthreads(num_threads);

    boost::scoped_array<boost::thread> threads(new boost::thread[num_threads]);
    for (int i = 0; i < num_threads; ++i) {
        boost::thread thr(boost::ref(D), i);
        threads[i].swap(thr);
    }

    for (int i = 0; i < num_threads; ++i)
        threads[i].join();
    
    // Back to one thread
    trans.set_nthreads(1);
    
}


// function to generate a simple graph
//
void generate_graph(Graph& g)
{

   // where is add_edge defined ?
   add_edge(0, 1, WeightedEdge(10), g);
   add_edge(0, 2, WeightedEdge(5), g);
   add_edge(1, 3, WeightedEdge(1), g);
   add_edge(1, 2, WeightedEdge(2), g);
   add_edge(2, 1, WeightedEdge(3), g);
   add_edge(2, 3, WeightedEdge(9), g);
   add_edge(2, 4, WeightedEdge(2), g);
   add_edge(3, 4, WeightedEdge(4), g);
   add_edge(4, 3, WeightedEdge(6), g);
   add_edge(4, 0, WeightedEdge(7), g);

}

void print_graph(Graph& g)
{
    graph_traits < Graph >::vertex_iterator i, end;
    graph_traits < Graph >::adjacency_iterator ai, a_end;
    property_map < Graph, vertex_index_t >::type
    index_map = get(vertex_index, g); // index_map, vertex_index  where these are defined ? how do we find there definitions ? how do we know those are there for us to use ?

    for (boost::tie(i, end) = vertices(g); i != end; ++i) {
        std::cout << get(index_map, *i);
        boost::tie(ai, a_end) = adjacent_vertices(*i, g);

        if (ai == a_end)
            std::cout << " has no children";
        else
            std::cout << " is the parent of ";
        for (; ai != a_end; ++ai) {
            std::cout << get(index_map, *ai);
            if (boost::next(ai) != a_end)
                std::cout << ", ";
        }
        std::cout << std::endl;
    }

    property_map<Graph, weight_type WeightedEdge::*>::type weights = get(&WeightedEdge::weight, g);

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
          std::cout << " --" << 1 << "-->" <<  get(weights, e) << "  ";
      }

      std::cout << std::endl;
    }


    /*graph_traits<Graph>::vertex_iterator iv, endv;
    graph_traits<Graph>::out_edge_iterator ei, edge_end;
    property_map<Graph, weight_type WeightedEdge::*>::type weights = get(&WeightedEdge::weight, g); // TODO whats going on here ?

    for (boost::tie(iv, endv) = vertices(g); iv != endv; ++iv) {
        std::cout << *iv << " ";
        for (boost::tie(ei, edge_end) = out_edges(*iv, g); ei != edge_end; ++ei)
            std::cout << " --" << *ei << "--> " << weights[*ei] << "  ";

        std::cout << std::endl;
    }*/
    std::cout << std::endl;
}

int test_main(int argc, char* argv[])
{

  std :: cout << "Initializing AM++ ...." << std :: endl;

  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true);
  amplusplus::transport trans = env.create_transport();

  std :: cout << "The first test starting ...." << std :: endl;

  Graph g(trans);
  generate_graph(g);
  //print_graph(g);
    

  return 0;
}
