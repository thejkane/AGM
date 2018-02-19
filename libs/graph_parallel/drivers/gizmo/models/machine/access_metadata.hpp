#ifndef __AGM_ACCESS_METADATA__
#define __AGM_ACCESS_METADATA__

class access_metadata {

public:
  uint64_t state_reads;
  uint64_t state_writes;
  uint64_t graph_data_reads;
  uint64_t wi_gens;

  access_metadata():state_reads(0),
		    state_writes(0),
		    graph_data_reads(0),
		    wi_gens(0){}

  void reset() {
    state_reads = 0;
    state_writes = 0;
    graph_data_reads = 0;
    wi_gens = 0;
  }

  inline uint64_t get_total_cost() const {
    // writes are treated as expensive than reads
    return (state_reads + graph_data_reads + (2*state_writes)
	    + (2*wi_gens));
  }
};

#endif
