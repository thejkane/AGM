#ifndef __AGM_LOGP_MODEL__
#define __AGM_LOGP_MODEL__

#include "../../utils/file_reader.hpp"

class distributed_logp {

private:
  int ranks;
  double latency;
  double overhead;
  double gap;

public:
  template<typename work_item_t>
  void initialize(transport& trans, const gizmo_config& configs) {
    double bwdth = configs.get_double(KEY_BANDWIDTH);
    assert(bwdth != 0);
    // Bandwidth is in MB/sec, convert it to bytes / sec
    gap = 1/bwdth;
    gap = gap / std::pow(2, 20);
    gap = gap * std::pow(10, 6); // gap in micro-seconds
    // gap is for 1 bytes
    // shouldnt we adjust value for n bytes ? TODO check


    latency = configs.get_double(KEY_LATENCY);
    ranks = configs.get_int(KEY_RANKS);

  }

};




#endif
