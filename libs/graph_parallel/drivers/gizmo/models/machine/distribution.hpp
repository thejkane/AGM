#ifndef __AGM_DISTRIBUTION__
#define __AGM_DISTRIBUTION__

enum graph_distribution_t {
  g_block
};


class block_graph_distribution {

private:
  uint64_t block_size;
  uint64_t num_vertices;
  int ranks;

public:
  block_graph_distribution() : num_vertices(0),
			       ranks(0),
			       block_size(0){
  }

  void initialize(uint64_t _n,
		  int _r) {

    num_vertices = _n;
    ranks = _r;
    block_size = (num_vertices / ranks) + ((num_vertices % ranks) > 0? 1 : 0);
  }

  template<typename Vertex>
  int owner(Vertex v) {
    assert(block_size != 0);
    return (v / block_size);
  }

};

#endif
