#ifndef __AGM_PRAM_MODEL__
#define __AGM_PRAM_MODEL__

#include "access_metadata.hpp"

enum operation_mode {
  pram_crcw,
  pram_ercw,
  pram_erew,
  pram_crew
};

class shared_memory {

private:
  int parallel_threads;
  int current_thread;
  uint64_t running_cost;
  operation_mode mode;
  uint64_t total_cost;

public:
  // programming model is always erew (cos the simulation is sequential)
  shared_memory(int threads):parallel_threads(threads), current_thread(0), running_cost(0),
	 mode(pram_erew), total_cost(0){}

  void update(const access_metadata& meta) {
    if (running_cost < meta.get_total_cost()) {
      running_cost = meta.get_total_cost();
    }

    ++current_thread;
    if (current_thread == parallel_threads) {
      total_cost += running_cost;
      running_cost = 0;
    }
  }

  void print() {
    std::cout << "[Shared-Memory] Total cost : " << total_cost << std::endl;
  }

};

#endif
