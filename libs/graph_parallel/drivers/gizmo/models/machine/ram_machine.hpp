#ifndef __AGM_RAM_MACHINE__
#define __AGM_RAM_MACHINE__

#include <stdio.h>
#include <queue>
#include "../../utils/file_reader.hpp"
#include "distribution.hpp"
#include "abstract_machine.hpp"


template<typename WorkItem,
	 typename StrictWeakOrderingRelation,
	 typename ProcessingFunction,
	 typename GraphDistribution=block_graph_distribution>
class ram_machine : public abstract_machine<WorkItem, StrictWeakOrderingRelation, 
						  ProcessingFunction, GraphDistribution> {

public:
  typedef abstract_machine<WorkItem, StrictWeakOrderingRelation,
			   ProcessingFunction, GraphDistribution> abstract_machine_t;

private:
  bool self_send;
  int latency_value;
  int send_overhead;
  int recv_overhead;

public:
  ram_machine(uint64_t n,
	      ProcessingFunction& _pf) : abstract_machine_t(n, _pf),
					 self_send(true),
					 latency_value(0),
					 send_overhead(0),
					 recv_overhead(0){}
  
  bool initialize(const gizmo_config& configs) {
    
    abstract_machine_t::initialize(configs);

    self_send = configs.get_bool(KEY_SELF_SEND_ENABLED);
    send_overhead = configs.get_int(KEY_SEND_OVERHEAD);
    recv_overhead = configs.get_int(KEY_RECV_OVERHEAD);
    latency_value = configs.get_int(KEY_LATENCY_VALUE);
  }

  int get_send_overhead() {
    // TODO decide precisely
    if (self_send)
      return 0;
    else 
      return 1;
  }

  int get_transmission_time(int source, int destination) {
    self_send = false;//debug
    if (self_send)
      return 0;
    else
      return 1700;
  }
};

struct ram_machine_gen {
public:
  template <typename WorkItem, typename StrictWeakOrderingRelation, typename ProcessingFunction>
  struct inner {
    typedef ram_machine<WorkItem, StrictWeakOrderingRelation, ProcessingFunction> type;
  };
};

#endif
