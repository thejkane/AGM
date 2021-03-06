// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

//======== Triangle Counting Algortihm================//
//===========================================================//


#ifndef BOOST_GRAPH_TC
#define BOOST_GRAPH_TC

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/parallel/iteration_macros.hpp> // for all_reduce
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <algorithm> // for std::min, std::max
#include <boost/format.hpp>
#include <iostream>
#include <atomic>
#include "boost/tuple/tuple.hpp"
#include "thread_pq_def.hpp"
#include <boost/graph/distributed/owner_defs.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

//for profiling
#ifdef CRAYPAT
#include <pat_api.h>
#endif

//===================== NOTE ===================================
// This datatype should never be used but is needed to build the
// distributed property map
//===================== NOTE ===================================
namespace amplusplus {
template <> struct make_mpi_datatype<std::vector<unsigned long long> > {MPI_Datatype get() const {return MPI_INT;}};
}

namespace amplusplus {
template <> struct make_mpi_datatype<std::vector<uint64_t> > {MPI_Datatype get() const {return MPI_INT;}};
}

namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
 typename IdDistribution,
 typename MessageGenerator = 
   amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class triangle_counting {

  typedef triangle_counting<Graph, IdDistribution, MessageGenerator> 
  self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;

  typedef std::vector<Vertex> InternalContainerType;

  typedef std::vector<InternalContainerType> PredecessorsType;
  typedef std::vector<InternalContainerType> SuccessorsType;

  typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
  typedef boost::iterator_property_map<typename PredecessorsType::iterator, VertexIndexMap>  PredecessorsMap;
  typedef boost::iterator_property_map<typename SuccessorsType::iterator, VertexIndexMap>  SuccessorsMap;

public:
  // < Target (Successor), Current, Predecessor of source >
  typedef std::pair<Vertex, std::pair<Vertex, Vertex> > work_item_t;

  // Current, Predecessor
  typedef std::pair<Vertex, Vertex> intermediate_data_t;

  static Vertex targetv(work_item_t wi) {
    return wi.first;
  }

  //  static Vertex differencev(work_item_t wi) {
  //  return wi.second;
  //}

  static Vertex predecessorv(work_item_t wi) {
    return wi.second;
  }

  //static level_t levelv(work_item_t wi) {
  //  return wi.second.second;
  //}


  struct processing_function;

  struct minimum_pair_first
  {
    template<typename T>
    const T& operator()(const T& x, const T& y) const { return x < y ? x : y; }

    template<typename F>
    struct result {
      typedef typename boost::function_traits<F>::arg1_type type;
    };
  };


  typedef typename MessageGenerator::template call_result<work_item_t, 
							  processing_function, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;

public:
  triangle_counting(Graph& g,
		    amplusplus::transport &t,
		    const IdDistribution& idd,
		    int offset,
		    MessageGenerator message_gen = 
		    MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
      g(g), 
      transport(t), 
      nthreads(t.get_nthreads()),
      id_distribution(idd),
      owner(get(vertex_owner, g)), 
      core_offset(offset),
      vec_predecessors(NULL),
      vertex_predecessors(NULL),
      vec_successors(NULL),
      vertex_successors(NULL),
  //#ifdef PRINT_DEBUG
    messages(0),
    totcounts(0),
  //#endif
    relax_msg(message_gen, transport, owner_from_pair<OwnerMap, work_item_t>(owner),
	      amplusplus::idempotent_combination(minimum_pair_first()))

  {
    initialize();
  }

  //destructor
  ~triangle_counting() {
    delete vec_predecessors;
    delete vertex_predecessors;

    delete vec_successors;
    delete vertex_successors;

#ifdef TRIANGLE_ENUMERATE
    for (int tid=0; tid < nthreads; ++tid) {
      work_item_t* arr = all_triangles[tid];
      delete [] arr;
    }

    delete [] all_triangles;
#endif

    delete [] threaded_triangle_indexes;
  }

  void operator() (int tid) { 
    run(tid); 
  }

  void run(int tid = 0);

  time_type get_start_time() {
    return start_time;
  }

  /*  std::vector<work_item_t>& get_thread_local_triangles(int tid) {
      return threaded_triangles[tid];
      }*/


  time_type get_elapsed_time() {
    return (end_time - start_time);
  }


#ifdef TRIANGLE_ENUMERATE
  // must be executed in a single thread
  void get_local_triangles(std::vector<work_item_t>& out) {
    for(int tid=0; tid < nthreads; ++tid) {
      out.insert(out.end(), all_triangles[tid], (all_triangles[tid]+threaded_triangle_indexes[tid]));
    }
  }
#endif

  uint64_t get_local_triangle_counts() {
    uint64_t total = 0;
    for(int tid=0; tid < nthreads; ++tid) {
      total += threaded_triangle_indexes[tid];
    }

    return total;
  }

  void print_triangle_counts() {
    std::cout << "========== Printing triangle counts per each thread ==============" << std::endl;
    for(int tid=0; tid < nthreads; ++tid) {
      std::cout << "[TID=" << tid << "] -- " << threaded_triangle_indexes[tid] << std::endl; 
    }
  }


private:
  void initialize();

  void process(const work_item_t& data, int tid);

  //  void process(uint64_t* buff, int sz, int tid);
  void handle_queue(const int tid, amplusplus::transport::end_epoch_request& request);

  template<typename SizeType>
  inline SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }

  template<typename SizeType>
  inline SizeType local_id(SizeType k) {
    return g.distribution().local(k);
  }


  work_item_t construct_wi(Vertex s,
			   Vertex c,
			   Vertex p) {
    // typedef std::pair<Vertex, std::pair<diff_t, level_t> > work_item_t;
    work_item_t wi(s, std::make_pair(c, p));
    return wi;
  }

  inline void copy_wi(work_item_t& to, const work_item_t& from) {
    to.first = from.first;
    to.second.first = from.second.first;
    to.second.second = from.second.second;
  }


  bool get_triangle(int tid, 
		    typename boost::graph_traits<Graph>::vertex_iterator& ite,
		    uint64_t& iteration,
		    uint64_t& processed,
		    work_item_t& wi);

/*
  bool get_triangle(int tid, 
		    typename boost::graph_traits<Graph>::vertex_iterator& ite,
		    uint64_t& iteration,
		    uint64_t& processed,
		    uint64_t* buffer,
		    int& offset);
*/

  const int dummy_first_member_for_init_order;
  const Graph& g;
  amplusplus::transport& transport;
  const int nthreads;
  const IdDistribution& id_distribution;
  const OwnerMap& owner; 
  int core_offset;
  PredecessorsType* vec_predecessors;
  PredecessorsMap* vertex_predecessors;
  SuccessorsType* vec_successors;
  SuccessorsMap* vertex_successors;

  RelaxMessage relax_msg;
  shared_ptr<amplusplus::detail::barrier> t_bar;

  time_type start_time;
  time_type end_time;

#ifdef TRIANGLE_ENUMERATE
  work_item_t** all_triangles;
#endif
  uint64_t* threaded_triangle_indexes;

  //#ifdef PRINT_DEBUG
  std::atomic<std::uint64_t> messages; 
  std::atomic<std::uint64_t> totcounts; 
  //#endif

};

#define TC_PARAMS                                   \
      typename Graph, typename IdDistribution, typename MessageGenerator

#define TC_TYPE                                    \
      triangle_counting<Graph, IdDistribution, MessageGenerator>


#ifdef TRIANGLE_ENUMERATE
#define MAX_TRIANGLE_COUNT 100000000
#endif

template<TC_PARAMS>
void
TC_TYPE::initialize() {

  relax_msg.set_handler(processing_function(*this));

  // Initialize predecessors per each vertex
  // creates a vector (size=0) per each vertex
  vec_predecessors = new PredecessorsType(num_vertices(g), std::vector<Vertex>(0)); 
  vertex_predecessors = new PredecessorsMap(vec_predecessors->begin(), 
					    get(boost::vertex_index, g));
  // Initialize successors per each vertex
  vec_successors = new SuccessorsType(num_vertices(g), std::vector<Vertex>(0));
  vertex_successors = new SuccessorsMap(vec_successors->begin(), get(boost::vertex_index, g));

  //  threaded_triangles.resize(nthreads);
  threaded_triangle_indexes = new uint64_t[nthreads];

  for (int i=0; i < nthreads; ++i) {
    threaded_triangle_indexes[i] = 0;
  }

#ifdef TRIANGLE_ENUMERATE
  all_triangles = new work_item_t*[nthreads];
#endif
}


template<TC_PARAMS>
void
TC_TYPE::process(const work_item_t& data, int tid) {
  Vertex s = data.first;
  Vertex p = data.second.second;
  
  if (std::binary_search((*vertex_predecessors)[s].begin(),
			 (*vertex_predecessors)[s].end(), p)) {

#ifdef TRIANGLE_ENUMERATE
    copy_wi(all_triangles[tid][threaded_triangle_indexes[tid]], data);
#endif
    ++(threaded_triangle_indexes[tid]);
  }
}


/*
class counting_output_iterator {
public:
  counting_output_iterator() : count{0} {}
  void operator++() { }
  counting_output_iterator& operator*() { return *this; }
  template<typename T>
  void operator=(T) { count++; }
  size_t get_count() { return count; }
private:
  size_t count;
};


template<TC_PARAMS>
void
TC_TYPE::process(uint64_t* buff, int sz, int tid) {
  Vertex s = buff[0];
  //int totpreds = sz - 1;

  counting_output_iterator output_ite;


  output_ite = std::set_intersection(buff+1, buff+sz,
			(*vec_predecessors)[s].begin(),
			(*vec_predecessors)[s].end(),
			output_ite);

  (threaded_triangle_indexes[tid]) += output_ite.get_count();
}
*/



//====================================================================================
// Threaded get_triangle -- Each thread call this function.
// Every candidate triangle is processed in a separate thread.
// tid -- Thread id
// ite -- Iterator of the vertex being processed (initially begin of vertices(g))
// iteration -- The iteration being processed by tid
// processed -- how many candidate triangle had being processed for previous vertices
// wi -- WorkItem to be processed
// Return value -- If true, wi is populated, else wi is not populated.
//====================================================================================
template<TC_PARAMS>
bool
TC_TYPE::get_triangle(int tid, 
		      typename boost::graph_traits<Graph>::vertex_iterator& ite,
		      uint64_t& iteration,
		      uint64_t& processed,
		      work_item_t& wi) {

  uint64_t element = nthreads*iteration + tid;

  uint64_t pred_count = (*vertex_predecessors)[*ite].size();
  uint64_t succ_count = (*vertex_successors)[*ite].size();

  uint64_t allcombinations = pred_count * succ_count;

  if ((element+1) > (allcombinations + processed)) {
    ++ite;
    processed = (allcombinations+processed);
    return false;
  } else {
    assert((processed <= element) && (element < (allcombinations + processed)));
    // find position within the current vertex
    uint64_t pos = element - processed;
    auto pred_index = pos / succ_count;
    auto succ_index = pos % succ_count;

    assert(succ_index < (*vertex_successors)[*ite].size());
    wi.first = (*vertex_successors)[*ite][succ_index];
    wi.second.first = (*ite);

    assert(pred_index < (*vertex_predecessors)[*ite].size());
    wi.second.second = (*vertex_predecessors)[*ite][pred_index];
    
    return true;
  }
}

/*
template<TC_PARAMS>
bool
TC_TYPE::get_triangle(int tid, 
		      typename boost::graph_traits<Graph>::vertex_iterator& ite,
		      uint64_t& iteration,
		      uint64_t& processed,
		      uint64_t* buffer,
		      int& offset) {

  uint64_t element = nthreads*iteration + tid;

  uint64_t pred_count = (*vertex_predecessors)[*ite].size();
  uint64_t succ_count = (*vertex_successors)[*ite].size();

  uint64_t allcombinations = pred_count * succ_count;

  if ((element+1) > (allcombinations + processed)) {
    ++ite;
    processed = (allcombinations+processed);
    return false;
  } else {
    assert((processed <= element) && (element < (allcombinations + processed)));
    // find position within the current vertex
    uint64_t pos = element - processed;
    auto pred_index = pos % pred_count;
    auto succ_index = pos / pred_count;

    offset = 0;
    assert(succ_index < (*vertex_successors)[*ite].size());
    Vertex succ = (*vertex_successors)[*ite][succ_index];
    buffer[offset] = succ;
    ++offset;

    for (int i=0; i < TC_BLOCK_SZ; ++i) {
      pred_index = i + pred_index;
      if (pred_index < pred_count) {
	Vertex pred = (*vertex_predecessors)[*ite][pred_index];
	buffer[offset] = pred;
	++offset;
      } else
	break;
    }

    iteration += offset;
    
    return true;
  }
}
*/

template<TC_PARAMS>
void
TC_TYPE::run(int tid) {
  AMPLUSPLUS_WITH_THREAD_ID(tid) {

    int nthreads = transport.get_nthreads();
    if (0 == tid) {
      // Set the number of threads to the barrier
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    }

    { amplusplus::scoped_epoch epoch(transport); }

    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
    t_bar->wait();

    // if two processes are running on the same node, core_offset
    // is important to achieve thread affinity
    if (pin(tid+core_offset) != 0) {
      std::cerr << "[ERROR] Unable to pin current thread to "
		<< "core : " << tid << std::endl;
      assert(false);
    }

    // wait till all threads are pinned
    t_bar->wait();
    { amplusplus::scoped_epoch epoch(transport); }

    validate_thread_core_relation();

    t_bar->wait();

#ifdef TRIANGLE_ENUMERATE
    uint64_t max_per_thread = MAX_TRIANGLE_COUNT/nthreads;
    for (int tid=0; tid < nthreads; ++tid) {
      all_triangles[tid] = new work_item_t[max_per_thread];
    }
#endif
    t_bar->wait();

    // should come before begin epoch
    start_time = get_time();
    // Start the algorithm

#ifdef CRAYPAT
    if (PAT_region_begin ( 1, "tcrun" ) == PAT_API_FAIL) {
      std::cout << "PAT begin failed ! " << std::endl;
      assert(false);
    }
#endif

    // Parallel iterate over all the vertices and collect predecessors and successors
    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	Vertex u = target(e, g);
	if (logical_id(v) > logical_id(u)) { // predecessors
	  (*vertex_predecessors)[v].push_back(u);
	} else {
	  (*vertex_successors)[v].push_back(u);
	}
      }

      std::sort((*vertex_predecessors)[v].begin(), (*vertex_predecessors)[v].end());
      std::sort((*vertex_successors)[v].begin(), (*vertex_successors)[v].end());
    }

    // #ifdef PRINT_DEBUG
    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
      totcounts += (uint64_t)((*vertex_predecessors)[v].size() * (*vertex_successors)[v].size());
    }
    t_bar->wait();
    if (tid == 0) 
      std::cout << "total width : " << totcounts.load() << std::endl;
    // #endif


#ifdef PRINT_DEBUG
    // Twitter graph specific
    Vertex v = 1803884;
    //    typename boost::graph_traits<Graph>::degree_size_type maxd;
    assert(out_degree(v, g) == 3081112);
    std::cout << "predecessors : " << (*vertex_predecessors)[v].size() << std::endl;
    std::cout << "successors : " << (*vertex_successors)[v].size() << std::endl;
    uint64_t totalmessages = ((*vertex_predecessors)[v].size() * (*vertex_successors)[v].size());
    std::cout << "total messages generated : " << totalmessages << std::endl;
#endif

    // need to make sure all predecessors and successors
    t_bar->wait();

    if (tid == 0)
      std::cout << "loop starting time : " << get_time() << std::endl;

#ifdef PRINT_DEBUG
    uint64_t viterations = 0;
    Vertex currentprocessing;
#endif

    uint64_t iteration = 0;
    uint64_t processed = 0;
    work_item_t wi;
    std::pair<typename boost::graph_traits<Graph>::vertex_iterator,
    	      typename boost::graph_traits<Graph>::vertex_iterator> itepair = vertices(g);
    
    typename boost::graph_traits<Graph>::vertex_iterator startite = itepair.first;

    {
      amplusplus::scoped_epoch epoch(transport);
      while (startite != itepair.second) {
	if (get_triangle(tid, startite,
			 iteration, 
			 processed,
			 wi)) {

	  relax_msg.send(wi);
	  //increase iteration for current thread
	  ++iteration;
	} 
      }
    }

/*
      int offset = 0;
      uint64_t* buffer = new uint64_t[TC_BLOCK_SZ];
    {
      amplusplus::scoped_epoch epoch(transport);
      while (startite != itepair.second) {
	if (get_triangle(tid, startite,
			 iteration, 
			 processed,
			 buffer,
			 offset)) {
	  process(buffer, offset, tid);
	  //relax_msg.send(wi);
	} 
      }
    }

    delete[] buffer;

*/

#ifdef PRINT_DEBUG
    std::cout << "messages for vetex v : " << v << " : " << messages.load() << std::endl;
#endif

#ifdef CRAYPAT
    if (PAT_region_end(1) == PAT_API_FAIL) {
      std::cout << "PAT end failed ! " << std::endl;
      assert(false);
    }
#endif

    t_bar->wait();
    end_time = get_time();
  }
}


template<TC_PARAMS>
struct TC_TYPE::
processing_function {
  
  processing_function() : self(NULL) {}
  processing_function(triangle_counting& self) : self(&self) {}
  
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
#ifdef PRINT_DEBUG
    std::cout << "Handler called in tid : " << tid << std::endl;
#endif

    self->process(data, tid);
  }

protected:
  triangle_counting* self;
};

}}}
#endif
