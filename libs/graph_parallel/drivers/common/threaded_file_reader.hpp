//============================================================================
// Copyright (C) 2017 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine
//============================================================================
#ifndef THREADED_FILE_READER
#define THREADED_FILE_READER

#include <boost/thread.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <dirent.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <boost/graph/util/utils.hpp>

enum graph_format {
  fmt_zero,
  fmt_one
};

struct vertex_comparator {
  template<typename Edge>
  bool operator()(const Edge& e1, const Edge& e2) {
    if (e1.first < e2.first)
      return true;
    else if (e1.first == e2.first) {
      return e1.second < e2.second;
    } else
      return false;
  }
};


template<typename vertices_size_type,
	 typename edge_t>
class segment_reader {

  typedef std::vector<edge_t> edges_t;
  typedef std::vector<edges_t> thread_edges_t;

public:
  segment_reader(std::vector<std::string>& sgmnts,
                 thread_edges_t& tedges,
                 int threads,
                 vertices_size_type st,
                 vertices_size_type cnt,
                 bool flp,
                 bool prep,
                 graph_format fmt):
    segments(sgmnts),
    threaded_edges(tedges),
    nthreads(threads),
    start(st),
    count(cnt),
    flip(flp),
    preprocess(prep),
    format(fmt)
  {}

  bool read_segment(std::string& fsegment,
                    int tid) {
    std::ifstream infile(fsegment.c_str());
    std::string line;
    while (std::getline(infile, line)) {
      if (line.compare(0, 1, "%") != 0) {
	vertices_size_type source, destination;
	std::istringstream iss(line);
	if ((iss >> source >> destination)) {
	  // do not read self-edges
	  if (preprocess && (source == destination))
	    continue;

	  if (format == fmt_one) {
	    assert(source > 0);
	    assert(destination > 0);

	    source = source - 1;
	    destination = destination - 1;
	  }

	  if (source >= start &&
	      (source < (start+count))) {
	    // source belongs to current locality
	    auto e = std::make_pair(source, destination); 
	    threaded_edges[tid].push_back(e);
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
	      threaded_edges[tid].push_back(e);
#ifdef PRINT_DEBUG
	      std::cout << "{(" << e.first << ", " << e.second << "), " << std::endl;
#endif
	    }
	  }
	} else {
	  std::cout << "[ERROR] An error occurred while reading graph"
		    << std::endl;
          infile.close();
	  return false;
	}
      }
    }

    infile.close();
    return true;
  }


  void operator()(int tid) {
    // each file is read in a different thread
    for (int i=tid; i < segments.size(); i=i+nthreads) {
#ifdef PRINT_DEBUG
      fprintf(stderr, "Thread %d, reading %s segment", tid, segments[i].c_str());
#endif
      read_segment(segments[i], tid);

      if (preprocess) {
        std::sort(threaded_edges[tid].begin(), 
                  threaded_edges[tid].end(), 
                  vertex_comparator());

        threaded_edges[tid].erase(std::unique(threaded_edges[tid].begin(), 
                                              threaded_edges[tid].end(), 
                                              [](edge_t e1, edge_t e2) { 
                                                return ((e1.first == e2.first) && 
                                                        (e1.second == e2.second));
                                              }), threaded_edges[tid].end());

      }
    }
  }

private:
  std::vector<std::string>& segments;
  thread_edges_t& threaded_edges;
  int nthreads;
  vertices_size_type start;
  vertices_size_type count;
  bool flip;
  bool preprocess;
  graph_format format;
};

template<typename vertices_size_type,
	 typename edge_t>
class threaded_file_reader {

  typedef std::vector<edge_t> edges_t;
  typedef std::vector<edges_t> thread_edges_t;

public:
  threaded_file_reader(const char* f, 
                       const char* segname,
                       int th) : fname(f),
                                 nthreads(th),
                                 segment_prefix(segname)
  {}

  threaded_file_reader(const char* f) : fname(f),
                                        nthreads(16)
  {}

  bool read(edges_t& alledges,
	    bool flip,
	    bool preprocess,
	    vertices_size_type& istart,
	    vertices_size_type& icount,
	    graph_format format);


private:
  std::string fname;
  int nthreads;
  std::string segment_prefix;

};

template<typename vertices_size_type,
	 typename edge_t>
bool threaded_file_reader<vertices_size_type, edge_t>::read(edges_t& alledges,
						    bool flip,
						    bool preprocess,
						    vertices_size_type& istart,
						    vertices_size_type& icount,
						    graph_format format) {
  std::vector<std::string> segments;
  DIR *pd;
  struct dirent *dir;
  pd = opendir(fname.c_str());

  if (NULL != pd) {
    while ((dir = readdir(pd)) != NULL) {
      std::string segment(dir->d_name);
      if (segment.compare(0, segment_prefix.length(), segment_prefix) == 0) {
        segment.insert(0, "/");
        segment.insert(0, fname);
        segments.push_back(segment);
      }
    }
  } else {
    std::cout << "[ERROR] Cannot open graph directory : " << fname << std::endl;
    return false;
  }

  closedir(pd);

#ifdef PRINT_DEBUG
  assert(segments.size() == 45);
  auto print = [](const std::string& n) { std::cout << " " << n; };
  std::for_each(segments.begin(), segments.end(), print);
  std::cout << std::endl;
#endif


  thread_edges_t threaded_edges(nthreads);

  segment_reader<vertices_size_type, edge_t> 
    segreader(segments, threaded_edges, nthreads,
	      istart,
	      icount,
	      flip,
	      preprocess,
	      format);

  boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
  for (int i = 0; i < nthreads - 1; ++i) {
    boost::thread thr(boost::ref(segreader), i + 1);
    threads[i].swap(thr);
  }

  segreader(0);

  for (int i = 0; i < nthreads - 1; ++i)
    threads[i].join();

  if (_RANK == 0)
    std::cout << "[INFO] Threads done reading..." << std::endl;

  std::size_t totalsz = 0;
  for(int k=0; k < nthreads; ++k) {
    //    std::cout << "Thread : " << k << " size : " << threaded_edges[k].size() << std::endl;
    totalsz += threaded_edges[k].size();
  }

  alledges.reserve(totalsz);

  //merging
  if (preprocess) {
    auto dfirst = alledges.end();
    for(int k=0; k < nthreads; ++k) {
      dfirst = std::merge(alledges.begin(), alledges.end(),
                 threaded_edges[k].begin(), threaded_edges[k].end(),
                 dfirst);
    }

    alledges.erase(std::unique(alledges.begin(), alledges.end(), 
			    [](edge_t e1, edge_t e2) { 
			      return ((e1.first == e2.first) && (e1.second == e2.second));
			    }), alledges.end());

  } else {
    for(int k=0; k < nthreads; ++k) {
      alledges.insert(alledges.end(), threaded_edges[k].begin(),
                      threaded_edges[k].end());
    }
  }

  if (_RANK == 0)
    std::cout << "Total edges read : " << alledges.size() << std::endl;

#ifdef PRINT_DEBUG
  if (!preprocess) {
    if (flip)
      assert(alledges.size() == (112307385 * 2));
    else
      assert(alledges.size() == 112307385);
  }
#endif

  return true;
  
}


#endif
