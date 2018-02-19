//============================================================================
// Copyright (C) 2017 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine
//============================================================================
#ifndef PBGL2_GRAPH_LOADER
#define PBGL2_GRAPH_LOADER
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <boost/lexical_cast.hpp>

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <inttypes.h>
#include <cstdlib>
#include <math.h>

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>
#include <am++/counter_coalesced_message_type.hpp>

#include <boost/graph/distributed/time_calc.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/use_mpi.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
#include <boost/graph/distributed/selector.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>
#include <boost/accumulators/accumulators.hpp>

#include <cds/init.h>       // for cds::Initialize and cds::Terminate
#include <cds/gc/hp.h>      // for cds::HP (Hazard Pointer) SMR

#include <random>
#include <parallel/algorithm>
#include <algorithm> // for std::min, std::max
#include <functional>
#include "threaded_file_reader.hpp"


#define HEADER_NODE_STR "# Nodes:"
#define HEADER_EDGE_STR "Edges:"

typedef uint64_t readerv_t;

class graph_reader_params {

public:
  graph_reader_params() : file_name(""),
			  read_header(false),
			  n(0),
			  m(0),
			  read_weights(false),
			  gen_weights(false),
			  no_weights(false),
			  preprocess(false),
			  read_graph(false),
			  graph_type("unknown"),
			  format(fmt_one),
			  flip(true),
			  segmented(false),
			  segment_prefix(""),
			  reader_threads(44)
  {}
  

  std::string file_name;
  bool read_header;
  readerv_t n;
  readerv_t m;
  bool read_weights;
  bool gen_weights;
  bool no_weights;
  bool preprocess;
  bool read_graph;
  std::string graph_type;
  graph_format format;
  bool flip;
  bool segmented;
  std::string segment_prefix;
  int reader_threads;
  bool mmio = false;

  bool parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];

      if (arg == "--file") {
	file_name = boost::lexical_cast<std::string>( argv[i+1]);
      }

      if (arg == "--format") {
	std::string sformat = boost::lexical_cast<std::string>( argv[i+1]);
	if (sformat == "zero")
	  format = fmt_zero;
      }

      if (arg == "--read-graph")
	read_graph = true;

      if (arg == "--segmented-file")
	segmented = true;

      if (arg == "--segment-prefix") {
	segment_prefix = boost::lexical_cast<std::string>( argv[i+1]);
      }

      if (arg == "--segment-reader-threads") {
	reader_threads = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--read-header") {
	read_header = true;
      }

      if (arg == "--vertices") {
	n = boost::lexical_cast<readerv_t>( argv[i+1]);
      }
      
      if (arg == "--edges") {
	m = boost::lexical_cast<readerv_t>( argv[i+1]); 
      }

      if (arg == "--flip-edge") {
	flip = boost::lexical_cast<bool>( argv[i+1]); 
      }

      if (arg == "--read-weight")
	read_weights = true;

      if (arg == "--gen-weight")
	gen_weights = true;

      if (arg == "--no-weight")
	no_weights = true;

      if (arg == "--preprocess")
	preprocess = true;

      if (arg == "--graph-type") {
	graph_type = boost::lexical_cast<std::string>( argv[i+1]);
	if (graph_type == "mmio")
	  mmio = true;
      }
    }

    if (read_graph) {
      FILE *file = fopen(file_name.c_str(), "r");
      if (file == NULL) {
	std::cout << "[ERROR] File cannot be open -- " << file_name << std::endl;
	return false;
      } else
	fclose(file);
    }

      // TODO more checking

    return true;
  }

  std::string get_graph_type() {
    return graph_type;
  }

  std::string get_generator_specific_input() {
    std::stringstream ss;
    ss << ", VERTICES :" << n
       << ", EDGES :" << m
       << ", PREPROCESS : " << preprocess;

    return ss.str();
  }


  bool load_graph_file() {
    return read_graph;
  }

  void print() {
    assert(_RANK != -1);
    if (_RANK == 0) {
      std::cout << "============= Printing Graph Generation Parameters =============" << std::endl;
      std::cout << "Graph File = " << file_name << std::endl;
      std::cout << "Graph Type = " << graph_type << std::endl;
      std::cout << "Read graph = " << read_graph << std::endl;
      std::cout << "Read Header = " << read_header << std::endl;
      std::cout << "Vertices to read  = " << n << std::endl;
      std::cout << "Edges to read = " << m << std::endl;
      std::cout << "Read weights = " << read_weights << std::endl;
      std::cout << "Generate weights = " << gen_weights << std::endl;
      std::cout << "No weights = " << no_weights << std::endl;
      std::cout << "Preprocess = " << preprocess << std::endl;
      std::cout << "Segmented Graph file = " << segmented << std::endl;
      std::cout << "Segmenet prefix = " << segment_prefix << std::endl;
      std::cout << "Segment file reader threads = " << reader_threads << std::endl;
      std::cout << "================================================================" << std::endl;
    }
  }
};


class graph_reader {

  typedef boost::compressed_sparse_row_graph<boost::directedS, boost::no_property, WeightedEdge, 
					     boost::no_property, boost::distributedS<readerv_t> > Digraph;
  typedef boost::graph_traits<Digraph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Digraph>::edges_size_type edges_size_type;
  typedef boost::graph_traits<Digraph>::vertices_size_type vertices_size_type;
  typedef boost::property_map<Digraph, boost::vertex_index_t>::type VertexIndexMap;
  typedef boost::property_map<Digraph, boost::vertex_owner_t>::const_type OwnerMap;
  typedef int32_t edge_weight_t;


public:
  graph_reader(amplusplus::transport& t,graph_reader_params& par):transport(t),
					       params(par), g(NULL)
  {}

  void load_graph_w_weights() {
    // TODO -- fill as necessary
  }

  void load_graph_w_generated_weights() {
    // TODO -- fill as necessary
  }

  /*  struct vertex_comparator {
    template<typename Edge>
    bool operator()(const Edge& e1, const Edge& e2) {
      if (e1.first < e2.first)
	return true;
      else if (e1.first == e2.first) {
	return e1.second < e2.second;
      } else
	return false;
    }
    };*/


  Digraph* graph() {
    return g;
  }

  // Preprocess the graph
  // 1. Remove self loops
  // 2. Remove parallel edges
  template<typename edge_t>
  void preprocess(std::vector<edge_t>& edges) {
    std::sort(edges.begin(), edges.end(), vertex_comparator());
    edges.erase(std::unique(edges.begin(), edges.end(), 
			    [](edge_t e1, edge_t e2) { 
			      return ((e1.first == e2.first) && (e1.second == e2.second));
			    }), edges.end());

#ifdef PRINT_DEBUG
    for (auto i=0; i < edges.size(); ++i) {
      std::cout << "(" << edges[i].first << ", " << edges[i].second << ")" << std::endl;
    }
#endif
  }


  void load_graph_no_weights() {
    if (_RANK == 0)
      std::cout << "[INFO] Reading graph without weights ..." << std::endl;

    vertices_size_type n; 
    edges_size_type m;
    // read header if vertices is not specified
    if (params.n == 0 || params.mmio) {
      read_header(n, m);
      params.n = n;
      params.m = m;
    } else {
      n = params.n;
      m = params.m;
    }

    uint64_t total_alloc = m;
    if (params.flip)
      total_alloc = m * (uint64_t)2;

    /*if (transport.size()==1) {
      std::cout << "[INFO] Vertex shuffling is enabled since there is only one process." << std::endl;
      // for shuffling
      shufflevector.resize(n);
      for(vertices_size_type j=0; j < n; ++j)
	shufflevector[j] = j;

      unsigned seed = 12345;
      std::shuffle(shufflevector.begin(), shufflevector.end(), std::default_random_engine(seed));
    } else {
      std::cout << "[INFO] Vertex shuffling is disabled" << std::endl;
      }*/

    typedef std::pair<Vertex, Vertex> edge_t;
    typedef std::vector<edge_t> edges_t;
    edges_t edges;
    //    std::cout << "Total allocations - " << total_alloc << std::endl;
    if (!params.segmented) {
      edges.reserve(total_alloc);
    }

    auto bdist = boost::parallel::block<vertices_size_type>(transport, n);
    // starting position for current rank
    vertices_size_type start_index = bdist.start(transport.rank());
    // locally stored ids in current rank
    vertices_size_type local_count = bdist.block_size(n);
    boost::parallel::variant_distribution<vertices_size_type> distrib = bdist;

#ifdef PRINT_DEBUG
    std::cout << "start_index=" << start_index << std::endl;
    std::cout << "local_count=" << local_count << std::endl;
    if (transport.size() == 1) {
      assert(local_count == n);
    }
#endif

    time_type rt1 = get_time();
    if (params.segmented) {
      threaded_file_reader<vertices_size_type, edge_t> 
	tfr(params.file_name.c_str(), 
	    params.segment_prefix.c_str(), 
	    params.reader_threads);

      tfr.read(edges, params.flip, params.preprocess,
	       start_index,
	       local_count,
	       params.format);
	       
    } else {

      if (params.mmio) {
	if (!read_mmio_edges(start_index, local_count, edges, params.flip)) {
	  std::cout << "[ERROR] Error while reading the graph" << std::endl;
	  assert(false);
	}

      } else {
	if (!read_edges(start_index, local_count, edges, params.flip)) {
	  std::cout << "[ERROR] Error while reading the graph" << std::endl;
	  assert(false);
	}
      }

      if (params.preprocess) {
	if (_RANK==0)
	  std::cout << "[INFO] Preprocessing the graph -- removing parallel edges" << std::endl;
	time_type pt1 = get_time();
	preprocess<edge_t>(edges);
	time_type pt2 = get_time();
	
	if (_RANK==0)
	  std::cout << "[INFO] The preprocessing time : " << pt2-pt1 << std::endl;
      }

    }

    time_type rt2 = get_time();
    if (_RANK==0)
      std::cout << "[INFO] Graph reading time : " << (rt2-rt1);

#ifdef PRINT_DEBUG
    std::cout << "edges.size()=" << edges.size() << ", total alloc=" << total_alloc << std::endl;
    //assert(edges.size() == total_alloc);
    std::cout << "read all edges " << std::endl;
#endif


    time_type gt1 = get_time();
    g = new Digraph(boost::edges_are_unsorted_multi_pass, edges.begin(), edges.end(),
	     n, transport, distrib);
    time_type gt2 = get_time();

    if (_RANK==0)
      std::cout << "[INFO] Local graph creation time : " << (gt2-gt1)
		<< std::endl;
    // Clear edge array above
    edges.clear();

#ifdef PRINT_DEBUG       
    Vertex maxu;
    boost::graph_traits<Digraph>::degree_size_type maxd = 0; 
    BGL_FORALL_VERTICES_T(v, *g, Digraph) {
      if (boost::out_degree(v, *g) > maxd) {
	maxd = boost::out_degree(v, *g);
	maxu = v;
      }
    }

    if (_RANK==0)
      std::cout << "Maximum outdegree : " << maxd << " max outdegree vertex : " << maxu << std::endl;
#endif

  }

  // MUST call finalize release memory associated with the graph
  void finalize() {
    delete g;
  }

  bool verify_graph() {
    uint64_t numv = num_vertices(*g);
    uint64_t nume = num_edges(*g);
    //    std::cout << "rank : " << transport.rank() << " nume : " << nume << std::endl;
    uint64_t allvs = 0;
    uint64_t alles = 0;
    MPI_Reduce(&numv, &allvs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nume, &alles, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);

    if (transport.rank() == 0) {
      //      std::cout << "rank : " << transport.rank() << " allvs : " << allvs << std::endl;
      //      std::cout << "rank : " << transport.rank() << " alles : " << alles << std::endl;
      assert((allvs == params.n) && "Loaded vertices are different from the input vertices");

      if (params.flip)
	assert((alles == (2 * params.m)) && "Loaded edges are different from the input edges");
      else
	assert((alles == params.m) && "Loaded edges are different from the input edges");
    }
  }

  void read_graph() {
    if (params.no_weights)
      load_graph_no_weights();
    else if (params.read_weights)
      load_graph_w_weights();
    else if (params.gen_weights)
      load_graph_w_generated_weights();
    else {
      std::cout << "[ERROR] Invalid graph reading option. Must specify"
		<< " one of --read-weight, --gen-weight or --no-weight"
		<< std::endl;
      assert(false);
    }

    { amplusplus::scoped_epoch epoch(transport); }

    if (!params.preprocess)
      verify_graph();

    if (transport.rank() == 0)
      std::cout << "[INFO] Graph loaded ..." << std::endl;
    
  }


private:
  //==========================
  // read the header of the graph file if exists
  //==========================
  void read_header(vertices_size_type& n, edges_size_type& m) {

    std::ifstream infile(params.file_name.c_str());
    std::string line;
    std::size_t node_str_len = strlen(HEADER_NODE_STR);

    while (std::getline(infile, line)) {
      if (line.compare(0, node_str_len, "%") != 0) {
	std::istringstream iss(line);
	vertices_size_type snum, dnum;
	if ((iss >> snum >> dnum >> m)) {
	  if (snum > dnum)
	    n = snum;
	  else
	    n = dnum;

	  std::cout << "[INFO] Vertices : " << n << ", Edges : " << m << std::endl;
	  break;
	}
      }
    }

    infile.close();
  }


  //==========================
  // edges read with out weights
  //=========================
  template<typename edges_t>
  bool read_mmio_edges(vertices_size_type start,
		  vertices_size_type count,
		  edges_t& edges,
		  bool flip = true){

    std::ifstream infile(params.file_name.c_str());
    std::string line;
    bool firstline = true; // execlude first line
    while (std::getline(infile, line)) {
      if (line.compare(0, 1, "%") != 0) {

	if (firstline) {
	  firstline = false;
	  continue;
	}

	vertices_size_type source, destination;
	std::istringstream iss(line);
	if ((iss >> source >> destination)) {
	  // do not read self-edges
	  if (params.preprocess && 
	      (source == destination))
	    continue;

	  if (params.format == fmt_one) {
	    assert(source > 0);
	    assert(destination > 0);

	    source = source - 1;
	    destination = destination - 1;
	  }

	  if (source >= start &&
	      (source < (start+count))) {
	    // source belongs to current locality
	    auto e = std::make_pair(source, destination); 
	    edges.push_back(e);
#ifdef PRINT_DEBUG
	    std::cout << "{(" << e.first << ", " << e.second << "), " << std::endl;
#endif
	  }
	  // if flip true flip the edge and read
	  if (flip) {
	    if (destination >= start &&
		(destination < (start+count))) {
	      // source belongs to current locality
	      auto e = std::make_pair(destination, source); 
	      edges.push_back(e);
#ifdef PRINT_DEBUG
	      std::cout << "{(" << e.first << ", " << e.second << "), " << std::endl;
#endif
	    }
	  }
	} else {
	  std::cout << "[ERROR] An error occurred while reading graph"
		    << std::endl;
	  return false;
	}
      }
    }
    return true;
  }


  vertices_size_type scramble(vertices_size_type v) {
    /*if (transport.size() == 1)
      return shufflevector[v];
    else
    return v;*/
    return v;
  }

  //==========================
  // edges read with out weights
  //=========================
  template<typename edges_t>
  bool read_edges(vertices_size_type start,
		  vertices_size_type count,
		  edges_t& edges,
		  bool flip = true){

    std::ifstream infile(params.file_name.c_str());
    std::string line;
    while (std::getline(infile, line)) {
      if (line.compare(0, 1, "%") != 0) {
	vertices_size_type source, destination;
	std::istringstream iss(line);
	if ((iss >> source >> destination)) {
	  // do not read self-edges
	  if (source == destination)
	    continue;

	  if (params.format == fmt_one) {
	    assert(source > 0);
	    assert(destination > 0);

	    source = source - 1;
	    destination = destination - 1;
	  }

	  if (source >= start &&
	      (source < (start+count))) {
	    // source belongs to current locality
	    auto e = std::make_pair(scramble(source), scramble(destination)); 
	    edges.push_back(e);
#ifdef PRINT_DEBUG
	    std::cout << "{(" << e.first << ", " << e.second << "), " << std::endl;
#endif
	  }
	  // if flip true flip the edge and read
	  if (flip) {
	    if (destination >= start &&
		(destination < (start+count))) {
	      // source belongs to current locality
	      auto e = std::make_pair(scramble(destination), scramble(source)); 
	      edges.push_back(e);
#ifdef PRINT_DEBUG
	      std::cout << "{(" << e.first << ", " << e.second << "), " << std::endl;
#endif
	    }
	  }
	} else {
	  std::cout << "[ERROR] An error occurred while reading graph"
		    << std::endl;
	  return false;
	}
      }
    }
    return true;
  }



  //==========================
  // edges read with weights
  //=========================
  template<typename edge_with_weight_t>
  bool read_edges_with_weight(vertices_size_type start,
			      vertices_size_type count,
			      edge_with_weight_t& edges){
    std::ifstream infile(params.file_name.c_str());
    std::string line;
    while (std::getline(infile, line)) {
      if (line.compare(0, 1, "#") != 0) {
	edges_size_type source, destination;
	edge_weight_t weight;
	std::istringstream iss(line);
	if ((iss >> source >> destination >> weight)) {
	  if (source >= start &&
	      (source < (start+count))) {
	    // source belongs to current locality
	    auto e = 
	      std::make_pair(source, std::make_pair(destination, weight)); 
	    edges.push_back(e);
#ifdef PRINT_DEBUG
	    std::cout << "{(" << e.first << ", " << e.second.first << "), " << weight
		      << "}" << std::endl;
#endif
	  }
	} else {
	  std::cout << "[ERROR] An error occurred while reading graph"
		    << std::endl;
	  return false;
	}
      }
    }
    return true;
  }


  //====================================
  // edges read while generating weights
  //====================================
  template<typename edges_size_type,
	   typename edge_with_weight_t,
	   typename weight_iterator_t>
  bool read_edges_gen_weight(vertices_size_type start,
			     vertices_size_type count,
			     edge_with_weight_t& edges,
			     weight_iterator_t wi) {
    std::ifstream infile(params.file_name.c_str());
    std::string line;
    while (std::getline(infile, line)) {
      if (line.compare(0, 1, "#") != 0) {
	edges_size_type source, destination;
	std::istringstream iss(line);
	if ((iss >> source >> destination)) {
	  if (source >= start &&
	      (source < (start+count))) {
	    // source belongs to current locality
	    auto e = 
	      std::make_pair(source, std::make_pair(destination, (*wi))); 
	    edges.push_back(e);
	    ++wi;
	  }
	} else {
	  std::cout << "[ERROR] An error occurred while reading graph"
		    << std::endl;
	  return false;
	}
      }
    }
    return true;
  }


private:
  amplusplus::transport& transport;
  graph_reader_params& params;
  Digraph* g;
  std::vector<vertices_size_type> shufflevector;
};

#endif
