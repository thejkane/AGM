#include <iostream>
#include "synthetic_generator.hpp"
#include "parser.hpp"
#include "executor.hpp"
#include <boost/graph/util/vertex_permutations.hpp>
#include <boost/graph/util/work_stats.hpp>
#include "ampp_perf_counters.hpp"


enum sssp_agm {
  delta,
  kla,
  dijkstra,
  chaotic
};

// TODO we need to get the full EAGM configuration
enum sssp_eagm {
  global,
  node,
  numa,
  thread
};

const uint64_t undefined_source = std::numeric_limits<uint64_t>::max();

class sssp_instance_params {
public:
  sssp_agm agm;
  sssp_eagm eagm;
  id_distribution_t id_distribution;
  bool verify;
  uint64_t source;
  bool level_sync;
  int delta_val;
  int k_val;
  int visited_threshold;
  uint64_t seed64;
  boost::rand48 synch_gen;
  boost::uniform_int<uint64_t>* prand_vertex;
  agm_work_stats work_stats;
  
  sssp_instance_params(int threads,
		       sssp_agm& _agm,
		       sssp_eagm _eagm,
		       id_distribution_t& idd,
		       bool v,
		       uint64_t s,
		       bool lsync,
		       int d,
                       int k,
		       int vthold,
		       uint64_t sd):
    agm(_agm), 
    eagm(_eagm), 
    id_distribution(idd), 
    verify(v),
    source(s),
    level_sync(lsync),
    delta_val(d),
    k_val(k),
    visited_threshold(vthold),
    seed64(sd),
    prand_vertex(NULL),
    work_stats(threads) {

    synch_gen.seed(seed64);
  }

  ~sssp_instance_params() {
    delete prand_vertex;
    prand_vertex = NULL;
  }

  std::string get_agm() {
    if (agm == delta)
      return "Delta";
    else if (agm == kla)
      return "KLA";
    else if (agm == chaotic)
      return "Chaotic";
    else if (agm == dijkstra)
      return "Dijkstra";
    else
      return "Invalid Algorithm";
  }

  std::string get_eagm() {
    if (eagm == global)
      return "Global";
    else if (eagm == node)
      return "Node";
    else if (eagm == numa)
      return "NUMA";
    else if (eagm == thread)
      return "Thread";
    else
      return "Invalid EAGM";
  }

  void print() {
    std::cout << "Algorithm : [" << get_agm()
	      << ", " << get_eagm() << "]" << std::endl;

    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "source : " << source << std::endl;
    std::cout << "level-sync? : " << level_sync << std::endl;
    std::cout << "delta val : " << delta_val << std::endl;
    std::cout << "kla val : " << k_val << std::endl;
    std::cout << "visited threshold : " << visited_threshold << std::endl;
    std::cout << "seed : " << seed64 << std::endl;
    std::cout << "verify : " << verify << std::endl;
  }

  std::string get_algorithm() {
    std::string algorithm = "[";
    algorithm += get_agm();
    algorithm += ", ";
    algorithm += get_eagm();
    algorithm += "]";

    return algorithm;
  }

  uint64_t get_source(uint64_t n) {
    if (source == undefined_source) {
      if (prand_vertex == NULL) {
	prand_vertex = new boost::uniform_int<uint64_t>(0, n-1);
      }
      // generate a random source
      return (*prand_vertex)(synch_gen);
    } else {
      return source;
    }
  }

  void print_summary() {
    work_stats.print_summary();
  }

};


class sssp_params {

private:
  std::vector<sssp_instance_params> instance_params;
  std::vector<sssp_agm> agms;
  std::vector<sssp_eagm> eagms;
  id_distribution_t id_distribution = horizontal; // default
  bool verify = false;
  std::vector<uint64_t> sources;
  std::vector<int> deltas;
  std::vector<int> ks;
  bool level_sync = false;
  int visited_threshold = 100;
  uint64_t seed64 = 12345;
  int threads = -1;

public:
  bool parse(int argc, char* argv[]){
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--agm") {
	if (strcmp(argv[i+1],"delta") == 0) {
	  agms.push_back(delta);
	} else if (strcmp(argv[i+1],"kla") == 0) {
	  agms.push_back(kla);
	} else if (strcmp(argv[i+1],"chaotic") == 0) {
	  agms.push_back(chaotic);
        } else if (strcmp(argv[i+1],"dijkstra") == 0) {
	  agms.push_back(dijkstra);
	} else {
	  std::cout << "[ERROR] Invalid AGM. Available : delta, choatic, kla" << std::endl;
	  return false;
	}
      }

      if (arg == "--threads") {
	threads = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--visited-threashold") {
	visited_threshold = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--seed") {
	seed64 = boost::lexical_cast<uint64_t>( argv[i+1] );
      }

      if (arg == "--source") {
	sources = extract_params<uint64_t>( argv[i+1] );
      }

      if (arg == "--delta-val") {
	deltas = extract_params<int>( argv[i+1] );
      }

      if (arg == "--kla-val") {
	ks = extract_params<int>( argv[i+1] );
      }

      if (arg == "--level-sync") {
	level_sync = true;
      }

      if (arg == "--eagm") {
	if (strcmp(argv[i+1],"global") == 0) {
	  eagms.push_back(global);
	} else if (strcmp(argv[i+1],"node") == 0) {
	  eagms.push_back(node);
	} else if (strcmp(argv[i+1],"numa") == 0) {
	  eagms.push_back(numa);
     	} else if (strcmp(argv[i+1],"thread") == 0) {
	  eagms.push_back(thread);
	} else {
	  std::cout << "[ERROR] Invalid EAGM. Available : global, node, numa, thread" << std::endl;
	  return false;
	}
      }

      if (arg == "--id-distribution") {
	if (strcmp(argv[i+1],"vertical") == 0)
	  id_distribution = vertical;
	else if (strcmp(argv[i+1],"horizontal") == 0)
	  id_distribution = horizontal;
	else {
	  std::cout << "Invalid id distribution type. Available types are vertical and horizontal" 
		    << std::endl;
	  return false;
	}
      }

      if (arg == "--verify") {
	verify = true;
      }
    }
  }

  void print() {
    /*    BOOST_FOREACH(sssp_agm a, agms) {
    BOOST_FOREACH(sssp_eagm ea, eagms) {
      std::cout << "Algorithm -- " << alg.get_algorithm() << std::endl;
    }
    }*/

    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "verify : " << verify << std::endl;
  }

  const std::vector<sssp_instance_params>&
  get_instance_params() {
    if (instance_params.empty()) {

      if (sources.empty())
	sources.push_back(undefined_source);
  
      if (deltas.empty()) {
	deltas.push_back(10);
      }

      if (ks.empty()) {
	ks.push_back(1);
      }

      if (threads == -1) {
	std::cout << "Number of threads must be specified using --threads ..." << std::endl;
	assert(false);
      }

      // at least one agm must be specified
      assert(!agms.empty());
      BOOST_FOREACH(sssp_agm a, agms) {
	BOOST_FOREACH (uint64_t s, sources) {
	  BOOST_FOREACH(int d, deltas) {
            BOOST_FOREACH(int k, ks) {
              if (!eagms.empty()) {
                BOOST_FOREACH(sssp_eagm ea, eagms) {
                  instance_params.push_back(sssp_instance_params(threads,
								 a, 
								 ea, 
                                                                 id_distribution, 
                                                                 verify, 
                                                                 s,
                                                                 level_sync,
                                                                 d,
                                                                 k,
                                                                 visited_threshold,
                                                                 seed64));
                }
              } else {
                // eagms not specified -- assume chaotic
                instance_params.push_back(sssp_instance_params(threads,
							       a, 
							       global, 
                                                               id_distribution, 
                                                               verify, 
                                                               s,
                                                               level_sync,
                                                               d,
                                                               k,
                                                               visited_threshold,
                                                               seed64));
              }
            }
          }
        }
      }
    }

    return instance_params;
  }
  
};
