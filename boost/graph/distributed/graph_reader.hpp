#ifndef PBGL2_GRAPH_READER
#define PBGL2_GRAPH_READER
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <boost/lexical_cast.hpp>

#define HEADER_NODE_STR "# Nodes:"
#define HEADER_EDGE_STR "Edges:"

template <typename block_node_t>
class graph_reader {

public:
  graph_reader(const char* filename):file_name(filename),
			     num_edges(0),
			     num_vertices(0){
    // TODO
  }

  void read_header() {
    std::ifstream infile(file_name.c_str());
    std::string line;
    std::size_t node_str_len = strlen(HEADER_NODE_STR);
    while (std::getline(infile, line)) {
      if (line.compare(0, node_str_len, HEADER_NODE_STR) == 0) {
	// we found nodes
	// find the position that edges starts
	std::size_t found = line.find(HEADER_EDGE_STR);
	if (found!=std::string::npos) {
	  // substring and extract nodes
	  std::string nstr = line.substr(strlen(HEADER_NODE_STR)+1, (found - node_str_len - 1));
	  //	  std::cout << "extracted nodes : " << nstr << std::endl;
	  num_vertices = std::stoull(nstr);

	  // extract number of edges
	  std::size_t e_start_point = found+strlen(HEADER_EDGE_STR)+1;
	  std::string estr = line.substr(e_start_point, (line.length()-e_start_point));
	  //	  std::cout << "extracted edges : " << estr << std::endl;
	  num_edges = std::stoull(estr);
	}
	break;
      }

    }
  }


  // edges read with weights
  template<typename edges_size_type,
	   typename edge_weight_t,
	   typename edge_with_weight_t>
  bool read_edges_with_weight(block_node_t start,
			      block_node_t count,
			      edge_with_weight_t& edges){
    std::ifstream infile(file_name.c_str());
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

  // edges read while generating weights
  template<typename edges_size_type,
	   typename edge_with_weight_t,
	   typename weight_iterator_t>
  bool read_edges_wo_weight(block_node_t start,
				block_node_t count,
			    edge_with_weight_t& edges,
			    weight_iterator_t wi) {

    std::ifstream infile(file_name.c_str());
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

  block_node_t get_num_edges() {
    return num_edges;
  }

  block_node_t get_num_vertices() {
    return num_vertices;
  }

private:
  std::string file_name;
  block_node_t num_edges;
  block_node_t num_vertices;
};

#endif
