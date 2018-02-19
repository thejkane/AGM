#ifndef __AGM_DISTRIBUTED_REAL__
#define __AGM_DISTRIBUTED_REAL__

#include <stdio.h>
#include <queue>
#include "../../utils/file_reader.hpp"
#include "distribution.hpp"
#include "abstract_machine.hpp"


template<typename WorkItem,
	 typename StrictWeakOrderingRelation,
	 typename ProcessingFunction,
	 typename GraphDistribution=block_graph_distribution>
class distributed_real_machine : public abstract_machine<WorkItem, StrictWeakOrderingRelation, 
						  ProcessingFunction, GraphDistribution> {

public:
  typedef abstract_machine<WorkItem, StrictWeakOrderingRelation,
			   ProcessingFunction, GraphDistribution> abstract_machine_t;

private:
  std::vector< std::vector<uint64_t> > latency_profile;
  bool latency_profile_enabled;
  int latency_value;
  int send_overhead;
  int recv_overhead;
  // MB/Sec
  int bandwidth;

public:
  distributed_real_machine(uint64_t n,
			   ProcessingFunction& _pf) : abstract_machine_t(n, _pf),
						      latency_profile_enabled(false),
						      latency_value(0),
						      send_overhead(0),
						      recv_overhead(0),
						      bandwidth(0){}
  
  bool initialize(const gizmo_config& configs) {
    
    abstract_machine_t::initialize(configs);

    latency_profile.resize(abstract_machine_t::get_ranks());
    for (int i=0; i < abstract_machine_t::get_ranks(); ++i) {
      latency_profile[i].resize(abstract_machine_t::get_ranks());
    }

    send_overhead = configs.get_int(KEY_SEND_OVERHEAD);
    recv_overhead = configs.get_int(KEY_RECV_OVERHEAD);
    bandwidth = configs.get_int(KEY_BANDWIDTH);
    latency_profile_enabled = configs.get_bool(KEY_LATENCY_PROFILE_ENABLED);
    
    if (latency_profile_enabled) {
      std::string file_name = configs.get(KEY_LATENCY_PROFILE);
      FILE *file = fopen(file_name.c_str(), "r");
      if (file == NULL) {
	std::cout << "[ERROR] File cannot be open -- " << file_name << std::endl;
	return false;
      } else
	fclose(file);

      std::cout << "[INFO] Reading latency profile from file : " << file_name << std::endl;
      std::ifstream infile(file_name.c_str());
      std::string line;

      while (std::getline(infile, line)) {
	std::istringstream iss(line);
	int i, j;
	uint64_t val;
	if ((iss >> i >> j >> val)) {
	  latency_profile[i][j] = val;
	}
      }

      infile.close();
    } else {
      latency_value = configs.get_int(KEY_LATENCY_VALUE);
    }
  }

  int get_latency(int src, int target) {
    if (latency_profile_enabled)
      return latency_profile[src][target];
    else 
      return latency_value;
  }

  int get_send_overhead() {
    return send_overhead;
  }

  int get_recv_overhead() {
    return recv_overhead;
  }

  int get_gap() {
    auto gap = (1 >> 20);
    gap = (gap / bandwidth) * sizeof(WorkItem);
    return gap;
  }

  int get_transmission_time(int source, int destination) {
    return (get_latency(source, destination) + get_recv_overhead() + 
	    get_gap());
  }

  void print() {
    std::cout << "[INFO] Printing latency profile ..." << std::endl;
    for (int i=0; i < abstract_machine_t::get_ranks(); ++i) {
      for (int j=0; j < abstract_machine_t::get_ranks(); ++j) {
	std::cout << i << "\t" << j << "\t" << latency_profile[i][j] << std::endl;
      }
    }
  }
};

struct distributed_real_machine_gen {
public:
  template <typename WorkItem, typename StrictWeakOrderingRelation, typename ProcessingFunction>
  struct inner {
    typedef distributed_real_machine<WorkItem, StrictWeakOrderingRelation, ProcessingFunction> type;
  };
};

#endif
