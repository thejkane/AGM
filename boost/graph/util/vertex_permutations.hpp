#ifndef PBGL_VERTEX_PERMUTATIONS
#define PBGL_VERTEX_PERMUTATIONS

#include <iostream>

//=================================================//
// vertex id distributions
enum id_distribution_t {
  vertical,
  horizontal,
  degree
};
//=================================================//


template<typename Graph>
struct graph_util {
  typedef typename boost::graph_traits<Graph>::vertices_size_type VerticesSzType;
public:
  static inline VerticesSzType local_id(const Graph g, VerticesSzType v) {
    return g.distribution().local(v);
  }
};


/**
* One-to-one mapping of vertices. In this case
* vertex is mapped to running id. Id is increase column wise.
* 1  4
* 2  5
* 3  6
*/
template<typename Graph>
struct block_id_distribution {
  typedef typename boost::graph_traits<Graph>::vertices_size_type VerticesSzType;
public:
  const static uint64_t local_id_mask = ((uint64_t)1 << 48)-1;

  block_id_distribution(const Graph& _pg, VerticesSzType& pn):g(_pg), n(pn) {
    if (_RANK == 0)
      std::cout << "[INFO] Using vertical id distribution" << std::endl;
  }

  
  template<typename SizeType>
  SizeType operator()(SizeType k) const {
    SizeType blocksz = g.distribution().block_size(0, n);
    return (get(get(boost::vertex_owner, g), k))*blocksz + g.distribution().local(k);
  }

  template<typename SizeType>
  amplusplus::transport::rank_type owner(SizeType v) const { // v is the logical id
    SizeType blocksz = g.distribution().block_size(0, n);
    return (v / blocksz);
  }

  template<typename SizeType>
  SizeType to_vertex_descriptor(SizeType v) const { // v is the logical id
    // TODO
    assert(false);
    return 0;
  }


private:
  const Graph& g;
  VerticesSzType n;
};


/**
* One-to-one mapping of vertices. In this case
* vertex is mapped to running id. Id is increase row wise.
* 1  2
* 3  4
* 5  6
*/
template<typename Graph>
struct row_id_distribution {
public:
  const static uint64_t local_id_mask = ((uint64_t)1 << 48)-1;

  row_id_distribution(const Graph& _pg, int ranks):g(_pg), totalranks(ranks) {
    if (_RANK == 0)
      std::cout << "[INFO] Using horizontal id distribution" << std::endl;
  }

  template<typename SizeType>
  SizeType operator()(SizeType k) const {

    SizeType localid = k & local_id_mask;
    SizeType owner = (k >> 48);

    auto offset = localid * totalranks;
    return offset + owner;
    //auto offset = g.distribution().local(k) * totalranks;
    //    return offset + get(get(boost::vertex_owner, g), k);
  }

  template<typename SizeType>
  amplusplus::transport::rank_type owner(SizeType v) const { // v is the logical id
    return (v % totalranks);
  }

  template<typename SizeType>
  SizeType to_vertex_descriptor(SizeType v) const { // v is the logical id
    SizeType localid = (v/totalranks);
    SizeType rank = (v % totalranks);
    
    return g.make_vertex_descriptor(rank, localid);
  }


private:
  const Graph& g;
  int totalranks;
};



/**
 * Degree based id distribution.
 */
template<typename Graph>
struct degree_id_distribution {
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
  typedef std::vector<uint64_t> VectorIndexMapType;
  typedef boost::iterator_property_map<VectorIndexMapType::iterator, VertexIndexMap> IndexMapType;


  struct degree_cmp {
  private:
    const Graph& g;
    bool ascending;

  public:
    degree_cmp(const Graph& pg, bool asc) : g(pg), ascending(asc) {}

    bool operator()(const Vertex& v1, const Vertex& v2) {
      if (ascending)
	return (boost::out_degree(v1, g) < boost::out_degree(v2, g));
      else
	return (boost::out_degree(v1, g) > boost::out_degree(v2, g));
    }
  };
  

public:
  degree_id_distribution(const Graph& _pg, int ranks, bool asc=true):g(_pg), totalranks(ranks) {
    if (_RANK == 0)
      std::cout << "[INFO] Using degree id distribution" << std::endl;

    auto numvs = num_vertices(g);
    ids.resize(numvs);
    indexmap = IndexMapType(ids.begin(), get(boost::vertex_index, g));

    VectorIndexMapType sortids;
    sortids.resize(numvs);
    Vertex sortind=0;

    if (_RANK == 0)
      std::cout << "[INFO] Sorting vertices by degree ..." << std::endl;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      sortids[sortind] = v;
      ++sortind;
    }

    std::sort(sortids.begin(), sortids.end(), degree_cmp(g, asc));

    /*sortind = 0;
    std::cout << "sortids=(";
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      std::cout << sortids[sortind] << "-" << boost::out_degree(sortids[sortind], g) << ", ";

      if (sortind != (num_vertices(g) - 1)) {
	assert(boost::out_degree(sortids[sortind], g) <= boost::out_degree(sortids[sortind+1], g));
      }

      ++sortind;
    }

    std::cout << ")" << std::endl;
    */

    if (_RANK == 0)
      std::cout << "[INFO] Calculating logical ids by degree ..." << std::endl;

    Vertex ordered = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      indexmap[sortids[ordered]] = ordered;
      ++ordered;
    }

    /*    ordered = 0;
    std::cout << "ordered=(";
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      std::cout << indexmap[v] << "-" << boost::out_degree(v, g) << ", ";
      ++ordered;
    }
    std::cout << ")" << std::endl;
    */
    
    sortids.clear();
    
    // verification


  }

  Vertex operator()(Vertex k) const {
    // TODO make it distributed
    return indexmap[k];
  }

  template<typename SizeType>
  SizeType to_vertex_descriptor(SizeType v) const { // v is the logical id
    // TODO
    assert(false);
    return 0;
  }

  template<typename SizeType>
  amplusplus::transport::rank_type owner(SizeType v) const { // v is the logical id
    return 0;
  }


private:
  const Graph& g;
  int totalranks;
  VectorIndexMapType ids;
  IndexMapType indexmap;
};



#endif
