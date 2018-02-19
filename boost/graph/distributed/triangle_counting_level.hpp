// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Triangle Counting Algortihm================//
//===========================================================//


#ifndef BOOST_GRAPH_TC_LEVEL
#define BOOST_GRAPH_TC_LEVEL

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/parallel/append_buffer.hpp>
#include <boost/parallel/buffer_iterator.hpp>
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

#define LEVEL_1 0
#define LEVEL_2 1
#define NEXT_LEVEL_1 2
#define TOTAL_BUCKETS 3

namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
	 typename FlagMap,
 typename IdDistribution,
 typename MessageGenerator = 
   amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class triangle_counting_level {

  typedef triangle_counting_level<Graph, IdDistribution, MessageGenerator> 
  self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;

  typedef append_buffer<Vertex, 10u> Bucket;
  typedef buffer_iterator<Vertex, 10u> BufferIterator;
  typedef typename std::vector<Bucket*>::size_type BucketIndex;

  typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

public:
  typedef uint32_t diff_t;
  typedef uint64_t VLevel;
  typedef int Level;

  static const uint64_t vertex_mask = std::numeric_limits<uint32_t>::max();

  // < Target , TriangleRoot, <Source, shouldaddtolevel1> >
  typedef std::pair<Vertex, Vertex> work_item_t;

  // intermediate work item_T
  typedef std::pair<Vertex, Vertex> intermediate_work_item_t;

  typedef append_buffer<intermediate_work_item_t, 10u> ActiveWorkSet;  

  inline static Vertex targetv(const work_item_t& wi) {
    return wi.first;
  }

  inline static Vertex sourcev(const work_item_t& wi) {
    return wi.second;
  }

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
  triangle_counting_level(Graph& g,
			  amplusplus::transport &t,
			  const IdDistribution& idd,
			  FlagMap flags,
			  int offset,			  
			  MessageGenerator message_gen = 
			  MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
      g(g), 
      transport(t), 
      nthreads(t.get_nthreads()),
      id_distribution(idd),
      owner(get(vertex_owner, g)), 
      level1flags(flags),
      core_offset(offset),
      level_1(LEVEL_1),
      next_level_1(NEXT_LEVEL_1),
    current_workset(0),
    next_workset(1),
    relax_msg(message_gen, transport, owner_from_pair<OwnerMap, work_item_t>(owner),
	      amplusplus::idempotent_combination(minimum_pair_first()))

  {
    initialize();
  }

  //destructor
  ~triangle_counting_level() {

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

  void print_data(Vertex vdiff) {
    diff_t v = vdiff >> 32;
    diff_t d = vdiff & vertex_mask;
    //    fprintf(stderr, "rank=%d,v=%zu,diff=%zu\n", transport.rank(), v, d);
  }

  time_type get_elapsed_time() {
    return (end_time - start_time);
  }


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

  inline static work_item_t construct_wi(Vertex target,
					 Vertex src) {

    return work_item_t(target, src);
  }


  inline void copy_wi(work_item_t& to, const work_item_t& from) {
    to.first = from.first;
    to.second = from.second;
  }


  bool get_triangle(int tid, 
		    typename boost::graph_traits<Graph>::vertex_iterator& ite,
		    uint64_t& iteration,
		    uint64_t& processed,
		    work_item_t& wi);


  const int dummy_first_member_for_init_order;
  const Graph& g;
  amplusplus::transport& transport;
  const int nthreads;
  const IdDistribution& id_distribution;
  const OwnerMap& owner; 
  FlagMap level1flags;
  int core_offset;
  std::vector<shared_ptr<Bucket> > buckets;
  std::vector<shared_ptr<ActiveWorkSet> > worksets;
  int level_1;
  int next_level_1;
  int current_workset;
  int next_workset;
  typename Bucket::size_type max_block_size = 20000000;

  RelaxMessage relax_msg;
  shared_ptr<amplusplus::detail::barrier> t_bar;

  time_type start_time;
  time_type end_time;

  uint64_t* threaded_triangle_indexes;

};

#define TC_PARAMS                                   \
      typename Graph, typename FlagMap, typename IdDistribution, typename MessageGenerator

#define TC_TYPE                                    \
      triangle_counting_level<Graph, FlagMap, IdDistribution, MessageGenerator>


template<TC_PARAMS>
void
TC_TYPE::initialize() {

  relax_msg.set_handler(processing_function(*this));

  // only three buckets
  // bucket 0 - root vertices of the DAG
  // bucket 1 - level-1 vertices of the DAG
  // bucket 2 - level-2 vertices of the DAG
  buckets.resize(TOTAL_BUCKETS);
  for (BucketIndex i = 0 ; i < buckets.size() ; ++i) {
    shared_ptr<Bucket> p(new Bucket);
    buckets[i].swap(p);
  }

  worksets.resize(2);
  for (auto i = 0 ; i < worksets.size() ; ++i) {
    shared_ptr<ActiveWorkSet> p(new ActiveWorkSet);
    worksets[i].swap(p);
  }

  BGL_FORALL_VERTICES_T(u, g, Graph) {
    // lower_neighbors are not initialized
    std::set<Vertex> locadjacencies;
    int count = 0;
    BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
      Vertex v = target(e, g);
      if (logical_id(v) < logical_id(u)) {
	if (locadjacencies.insert(v).second) { //to avoid parallel edges
	  ++count;
	}
      }
    }

    if (count == 0) {
      if (out_degree(u, g) > 1) {
	//	fprintf(stderr, "rank=%d,root=%zu", transport.rank(), u);
	intermediate_work_item_t iwi(u, u);
	worksets[current_workset]->push_back(iwi);
      }
    }

  }


  //  threaded_triangles.resize(nthreads);
  threaded_triangle_indexes = new uint64_t[nthreads];

  for (int i=0; i < nthreads; ++i) {
    threaded_triangle_indexes[i] = 0;
  }

}


template<TC_PARAMS>
void
TC_TYPE::process(const work_item_t& data, int tid) {
  using boost::parallel::bool_compare_and_swap;

  Vertex src = sourcev(data);
  Vertex currv = targetv(data);  

  uint64_t w = (uint32_t)src;
  w = w << 32;
  w = w | (uint32_t)currv;
  buckets[LEVEL_2]->push_back(w);
}

template<TC_PARAMS>
void
TC_TYPE::run(int tid) {

  using boost::parallel::bool_compare_and_swap;

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

    // should come before begin epoch
    start_time = get_time();
    // Start the algorithm

#ifdef CRAYPAT
    if (PAT_region_begin ( 1, "tcrun" ) == PAT_API_FAIL) {
      std::cout << "PAT begin failed ! " << std::endl;
      assert(false);
    }
#endif

    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	Vertex u = target(e, g);
	if (logical_id(v) > logical_id(u)) { // predecessors
	  uint64_t w = (uint32_t)u;
	  w = w << 32; 
	  w = w | (uint32_t)v;
	  buckets[level_1]->push_back(w);
	}
      }
    }

    time_type l1tme = get_time();
    if (tid == 0)
      std::cout << "Time to fill L1 buckets : " << (l1tme - start_time) << std::endl;

    t_bar->wait();

    time_type feptime = get_time();

    {
      amplusplus::scoped_epoch epoch(transport);
      typename Bucket::size_type current_bucket_size = buckets[level_1]->size();
      for (typename Bucket::size_type i = tid ; 
	   i < current_bucket_size ; i+= nthreads) {

	uint64_t w = (*buckets[level_1])[i];
	Vertex u = w >> 32;
	Vertex v = w & vertex_mask;
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex s = target(e, g);
	  if (logical_id(s) > logical_id(v)) { // predecessors
	    work_item_t wi = construct_wi(s, u);
	    relax_msg.send(wi);
	  }
	}
      }      
    }

    time_type efeptime = get_time();

    if (tid == 0)
      std::cout << "Time to fill L2 bucket : " << efeptime-feptime << std::endl;


    typename Bucket::size_type block_sz_1 = max_block_size;
    typename Bucket::size_type block_sz_2 = max_block_size;

    auto level1sz = buckets[level_1]->size();
    auto level2sz = buckets[LEVEL_2]->size();


    time_type tint1 = get_time();
    // do set intersection
    if ((level1sz != 0) &&
	(level2sz != 0)) {

      if ((level1sz/nthreads) < block_sz_1) {
	block_sz_1 = buckets[level_1]->size()/nthreads;

	if (block_sz_1 == 0)
	  block_sz_1 = 1;
      }

      if ((level2sz/nthreads) < block_sz_2) {
	block_sz_2 = buckets[LEVEL_2]->size()/nthreads;

	if (block_sz_2 == 0)
	  block_sz_2 = 1;
      }

      auto numblocks_1 = (level1sz + block_sz_1 - 1)/block_sz_1;
      auto numblocks_2 = (level2sz + block_sz_2 - 1)/block_sz_2;

      // sort second level blocks in parallel
      for (auto pos = tid; pos < numblocks_2; pos+=nthreads) {
	typename Bucket::size_type begin = pos * block_sz_2;
	typename Bucket::size_type size = std::min(block_sz_2, (level2sz-(pos*block_sz_2)));
	typename Bucket::size_type end = begin+size;

	BufferIterator level2bufbegin(buckets[LEVEL_2].get(), begin);
	BufferIterator level2bufend(buckets[LEVEL_2].get(), end);
	  
	std::sort(level2bufbegin, level2bufend);
      }


      t_bar->wait();

      for (auto pos = tid; pos < numblocks_1; pos+=nthreads) {
	typename Bucket::size_type begin = pos * block_sz_1;
	const typename Bucket::size_type fixedbegin = begin;
	typename Bucket::size_type size = std::min(block_sz_1, (level1sz-(pos*block_sz_1)));
	typename Bucket::size_type end = begin+size;
	//fprintf(stderr, "b=%zu,e=%zu,tid=%zu", begin, end, tid);

	// sort the level-1 buffer
	BufferIterator level1bufbegin(buckets[level_1].get(), begin);
	BufferIterator level1bufend(buckets[level_1].get(), end);
	std::sort(level1bufbegin, level1bufend);

	// Level 2 items
	for (auto ipos = 0; ipos < numblocks_2; ++ipos) {
	  typename Bucket::size_type ibegin = ipos * block_sz_2;
	  typename Bucket::size_type isize = std::min(block_sz_2, (level2sz-(ipos*block_sz_2)));
	  typename Bucket::size_type iend = ibegin+isize; 
	  //	    fprintf(stderr, "b2=%zu,e2=%zu,tid=%zu", ibegin, iend, tid);
	  // do the actual matchine
	  while((begin != end) && (ibegin != iend)) {
	    if ((*buckets[level_1])[begin] < (*buckets[LEVEL_2])[ibegin]) {
	      ++begin;
	    } else if ((*buckets[LEVEL_2])[ibegin] < (*buckets[level_1])[begin]) {
	      ++ibegin;
	    } else {
	      ++(threaded_triangle_indexes[tid]);

	      while(((ibegin+1) != iend) &&
		    ((*buckets[LEVEL_2])[ibegin] == (*buckets[LEVEL_2])[ibegin+1])) {
      		++(threaded_triangle_indexes[tid]);
		++ibegin;
	      }

	      ++begin;
	      ++ibegin;
	    }
	  }
	    
	  // End of processing one Level 2 bucket
	  // Now we need to reset the starting point level-1 bucket to process next level2 bucket
	  begin = fixedbegin;

	}
      }
    }

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
  processing_function(triangle_counting_level& self) : self(&self) {}
  
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
#ifdef PRINT_DEBUG
    std::cout << "Handler called in tid : " << tid << std::endl;
#endif

    self->process(data, tid);
  }

protected:
  triangle_counting_level* self;
};

}}}
#endif
