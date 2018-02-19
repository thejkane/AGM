#include "./common/synthetic_generator.hpp"
#include "./common/graph_loader.hpp"
#include "./common/parser.hpp"

#include "gizmo/utils/utils.hpp"
#include "gizmo/utils/file_reader.hpp"
#include "gizmo/stats/stat.hpp"
#include "gizmo/applications/bfs.hpp"
#include "gizmo/models/machine/shared_memory.hpp"
#include "gizmo/models/machine/distributed_real.hpp"
#include "gizmo/models/machine/ram_machine.hpp"


#include <limits.h>

void* util::p_callback = NULL;

struct swo_1 {
public:
  bool operator()(int i, int j) {
    return (i < j);
  }
};

//chaotic
struct swo_2 {
public:
  bool operator()(int i, int j) {
    return false;
  }
};

//similar to delta
struct swo_3 {
public:
  bool operator()(int i, int j) {
    return ((i/3) < (j/3));
  }
};

class test_buckets {

public:
  void run_simple_integrated_tests() {
    test_buckets tbuckets;
    tbuckets.run_tests();

    typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS,
				    boost::no_property, boost::property < boost::edge_weight_t, int > > graph_t;
    typedef std::pair<int, int> Edge;

    const int num_nodes = 5;
    Edge edge_array[] = { Edge(1, 3), Edge(2, 2), Edge(2, 4), Edge(2, 5),
			  Edge(3, 2), Edge(3, 4), Edge(4, 5), Edge(5, 1), Edge(5, 2)
    };
    int weights[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    int num_arcs = sizeof(edge_array) / sizeof(Edge);

    graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);

    stat_reader sr;
    gizmo_config gizmocfg;
    bfs_simulator< graph_t> bfssim;
    bfssim.template simulate<chaotic_ordering_gen, distributed_real_machine_gen>(g, sr, gizmocfg);
  }

  void run_test_1() {
    buckets<int, swo_1> all_buckets;
    for (int i=0; i < 10; ++i) {
      all_buckets.push(i);
    }

    buckets<int, swo_1>::bucket_iterator_t bbegin, bend;
    int j = 0;
    while(!all_buckets.empty()) {
      std::tie(bbegin, bend) = all_buckets.top_bucket();
      std::cout << "Bucket :" << j << " - {";
      while (bbegin != bend) {
        std::cout << *bbegin << ", ";
        ++bbegin;
      }
      std::cout << "}" << std::endl;
      all_buckets.pop_bucket();
      ++j;
    }
    
  }

  void run_test_2() {
    buckets<int, swo_2> all_buckets;
    for (int i=0; i < 10; ++i) {
      all_buckets.push(i);
    }

    buckets<int, swo_1>::bucket_iterator_t bbegin, bend;
    int j = 0;
    while(!all_buckets.empty()) {
      std::tie(bbegin, bend) = all_buckets.top_bucket();
      std::cout << "Bucket :" << j << " - {";
      while (bbegin != bend) {
        std::cout << *bbegin << ", ";
        ++bbegin;
      }
      std::cout << "}" << std::endl;
      all_buckets.pop_bucket();
      ++j;
    }

    //all_buckets.print_buckets();
  }

  void run_test_3() {
    buckets<int, swo_3> all_buckets;
    for (int i=0; i < 10; ++i) {
      all_buckets.push(i);
    }

    buckets<int, swo_1>::bucket_iterator_t bbegin, bend;
    int j = 0;
    while(!all_buckets.empty()) {
      std::tie(bbegin, bend) = all_buckets.top_bucket();
      std::cout << "Bucket :" << j << " - {";
      while (bbegin != bend) {
        std::cout << *bbegin << ", ";
        ++bbegin;
      }
      std::cout << "}" << std::endl;
      all_buckets.pop_bucket();
      ++j;
    }
  }


  void run_tests() {
    info("Running Tests : Every integer is a separate equivalence class");
    run_test_1();
    info("Running Tests : All integers are in the same equivalence classs (Chaotic)");
    run_test_2();
    info("Running Tests : All integers are in three equivalence classses (~Delta)");
    run_test_3();

  }
};


template<typename Graph>
void simulate(Graph& g, gizmo_config& config) {
  info("Invoking simulation ...");
  stat_reader sr;

  bfs_simulator<Graph> bfssim;

  time_type start = get_time();
  if (config.get_machine_model() == distributed_real) {
    info("Invoking simulation on distributed real machine model ...");
    bfssim.template simulate<chaotic_ordering_gen, 
			     distributed_real_machine_gen>(g, sr, config);
  } else if (config.get_machine_model() == ram) {
    info("Invoking simulation on ram machine model ...");
    bfssim.template simulate<chaotic_ordering_gen, 
			     ram_machine_gen>(g, sr, config);

  } else {
    std::cout << "[ERROR] Invalid machine model." << std::endl;
    return;
  }
  time_type end = get_time();

  std::cout << "[INFO] Simulation done in : "<< (end-start) << " sec. " << std::endl;

  sr.print();
}

uint64_t execution_time = 0;

int main(int argc, char* argv[]) {

  info("Starting AGM simulator ...");

  runtime_config_params rcp;
  rcp.parse(argc, argv);

  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true,
							    1/*recvdepth*/,
							    1/*polls*/,
							    rcp.get_flow_control());
  amplusplus::transport trans = env.create_transport();

  _RANK = trans.rank();

  rcp.print();

  gizmo_config gizmocfg;
  config_file_reader cfr;
  if (cfr.parse(argc, argv))
    cfr.read_configurations(gizmocfg);

  gizmocfg.print();

  // Are we reading a graph from a file ?
  graph_reader_params readeparams;
  if (!readeparams.parse(argc, argv))
    return -1;

  // No we are not -- then generate the graph
  if (!readeparams.read_graph) {
    synthetic_generator sg(trans);

    graph_gen_params gparams;
    if (!gparams.parse(argc, argv))
      return -1;

    gparams.print();
    sg.generate_graph(gparams);
      
    if(_RANK == 0)
      std::cout << "[INFO] Done generating the graph ..." << std::endl;

    simulate(*sg.graph(), gizmocfg);
    sg.finalize();

  } else {
    graph_reader gr(trans, readeparams);
    // we are reading a graph from a file
    readeparams.print();
    gr.read_graph();

    simulate(*gr.graph(), gizmocfg);
    gr.finalize();
  }
}
