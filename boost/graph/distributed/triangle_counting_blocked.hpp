// Copyright (C) 2015-2016 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Triangle Counting Algortihm================//
//===========================================================//


#ifndef BOOST_GRAPH_TC_BLOCKED
#define BOOST_GRAPH_TC_BLOCKED

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

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
#include <am++/size_coalesced_message_type.hpp>
#include <chrono>
#include <unordered_map>
#include <boost/functional/hash.hpp>

//for profiling
#ifdef CRAYPAT
#include <pat_api.h>
#endif


namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
	 typename IdDistribution, 
	 typename NeighborMap,         
	 typename MessageGenerator = 
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> > 
class triangle_counting_blocked {

  typedef triangle_counting_blocked<Graph, IdDistribution, NeighborMap, MessageGenerator> self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;

  typedef typename std::pair< typename std::vector<Vertex>::iterator,
			      typename std::vector<Vertex>::iterator > IteratorPair_t;

  struct block_msg;
  struct processing_function;
  // AM++ message type
  //  typedef amplusplus::message_type<Vertex> tm_type;
  typedef typename amplusplus::size_coalesced_message_type_gen::inner<block_msg, processing_function>::type tm_type;

  struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
      return boost::hash_value(p);
    }
  };

  typedef std::unordered_map<std::pair<Vertex, Vertex>, bool, pair_hash> HashMap_t;


public:

  struct degree_processing_function;
  typedef std::pair<Vertex, std::pair<Vertex, uint32_t > > work_item_t;
  struct minimum_pair_first
  {
    template<typename T>
    const T& operator()(const T& x, const T& y) const { return x.first < y.first ? x : y; }

    template<typename F>
    struct result {
      typedef typename boost::function_traits<F>::arg1_type type;
    };
  };


  typedef typename MessageGenerator::template call_result<work_item_t, 
							  degree_processing_function, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;



  static Vertex targetv(work_item_t wi) {
    return wi.first;
  }

public:
  triangle_counting_blocked(Graph& g,
			    amplusplus::transport &t,
			    const IdDistribution& idd,
			    int offset,
			    uint64_t& bz,
			    uint64_t& sucbz,
			    uint64_t& csz,
			    NeighborMap& pred,
			    NeighborMap& succ,
			   MessageGenerator message_gen =
			    MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 17)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<work_item_t>(), 0)),
      g(g), 
      transport(t), 
      nthreads(t.get_nthreads()),
      id_distribution(idd),
      owner(get(vertex_owner, g)), 
      core_offset(offset),
      block_size(bz),
      suc_block_size(sucbz),
    coalescing_size(csz),
    vertex_predecessors(pred),
    vertex_successors(succ),
    msg_type(amplusplus::size_coalesced_message_type_gen(coalescing_size), t),
    relax_msg(message_gen, transport, owner_from_pair<OwnerMap, work_item_t>(owner),
	      amplusplus::idempotent_combination(minimum_pair_first()))

  {
    initialize();
  }

  //destructor
  ~triangle_counting_blocked() {
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
      std::cout << "[Rank=" << transport.rank() << "TID=" << tid << "] -- " << threaded_triangle_indexes[tid] << std::endl; 
    }
  }

  void write_to_file() {
    std::string sfile = "succ-prof-rank-";
    std::string r = std::to_string((int)transport.rank());
    sfile += r;
    sfile += ".txt";

    std::ofstream outfile;
    outfile.open(sfile.c_str());
    outfile << "Vertex" << "\t" << "Successors" << std::endl; 

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      outfile << logical_id(v) << "\t" << vertex_successors[v].size() << std::endl;      
    }
    
    outfile.close();

    sfile = "pred-prof-rank-";
    sfile += r;
    sfile += ".txt";

    outfile.open(sfile.c_str());
    outfile << "Vertex" << "\t" << "Predecessors" << std::endl; 

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      outfile << logical_id(v) << "\t" << vertex_predecessors[v].size() << std::endl;      
    }
    
    outfile.close();

  }

#ifdef TC_STATS
  void print_stats() {
    uint64_t all_send_msgs = 0;
    uint64_t all_recv_msgs = 0;
    uint64_t all_local_msgs = 0;
    uint64_t all_preds = 0;
    uint64_t all_succs = 0;
    uint64_t all_succ_preds = 0;
    uint64_t tot_avg_set_int_times = 0;
    uint64_t total_set_inters = 0;
    uint64_t all_saved_msgs = 0;
    uint64_t all_comparisons = 0;
    uint64_t all_bytes_over_nw = 0;
    uint64_t all_init_succs = 0;
    uint64_t all_init_preds = 0;
    uint64_t all_predicted_psp_bytes = 0;
    uint64_t all_predicted_sps_bytes = 0;
    uint64_t vmax_pred_block = 0;
    uint64_t vmax_succ_block = 0;
    

    uint64_t t_all_send_msgs = 0;
    uint64_t t_all_recv_msgs = 0;
    uint64_t t_all_local_msgs = 0;
    uint64_t t_all_preds = 0;
    uint64_t t_all_succs = 0;
    uint64_t t_all_succ_preds = 0;
    uint64_t t_tot_avg_set_int_times = 0;
    uint64_t t_total_set_inters = 0;
    uint64_t t_all_saved_msgs = 0;
    uint64_t t_all_comparisons = 0;
    uint64_t t_all_bytes_over_nw = 0;
    uint64_t t_all_init_succs = 0;
    uint64_t t_all_init_preds = 0;
    uint64_t t_all_predicted_psp_bytes = 0;
    uint64_t t_all_predicted_sps_bytes = 0;
    uint64_t t_vmax_pred_block = 0;
    uint64_t t_vmax_succ_block = 0;


    for (int i=0; i < nthreads; ++i) {
      t_all_send_msgs += send_msgs[i];
      t_all_recv_msgs += recv_msgs[i];
      t_all_local_msgs += local_msgs[i];
      t_all_preds += total_preds[i];
      t_all_succs += total_succs[i];
      t_all_succ_preds += total_succs_preds[i];
      t_tot_avg_set_int_times += tot_set_int_time[i];
      t_total_set_inters += total_setints[i];
      t_all_saved_msgs += saved_msgs[i];
      t_all_comparisons += tot_comparisons[i];
      t_all_bytes_over_nw += num_bytes_sent_over_nw[i];
      t_all_init_succs += init_succs[i];
      t_all_init_preds += init_preds[i];
      t_all_predicted_psp_bytes += predicted_psp_bytes[i];
      t_all_predicted_sps_bytes += predicted_sps_bytes[i];

      if (max_pred_block[i] > t_vmax_pred_block)
	t_vmax_pred_block = max_pred_block[i];

      if (max_succ_block[i] > t_vmax_succ_block)
	t_vmax_succ_block = max_succ_block[i];

    }

    MPI_Reduce(&t_all_send_msgs, &all_send_msgs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_recv_msgs, &all_recv_msgs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_local_msgs, &all_local_msgs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_preds, &all_preds, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_succs, &all_succs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_succ_preds, &all_succ_preds, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_tot_avg_set_int_times, &tot_avg_set_int_times, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_total_set_inters, &total_set_inters, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_saved_msgs, &all_saved_msgs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_comparisons, &all_comparisons, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_bytes_over_nw, &all_bytes_over_nw, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_init_succs, &all_init_succs, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_init_preds, &all_init_preds, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_predicted_psp_bytes, &all_predicted_psp_bytes, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_all_predicted_sps_bytes, &all_predicted_sps_bytes, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_vmax_pred_block, &vmax_pred_block, 
	       1, MPI_LONG_LONG_INT , MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_vmax_succ_block, &vmax_succ_block, 
	       1, MPI_LONG_LONG_INT , MPI_MAX, 0, MPI_COMM_WORLD);


    auto nranks = transport.size();
    if (transport.rank() == 0) {
      std::string with_ordering = "Yes";
#ifdef TC_NO_ORDERING
      with_ordering = "No";
#endif
      std::cout << "[INFO][STATS] Per rank stats. Ranks : " << nranks
		<< ", Total pre-calculated successors : " << all_init_succs
		<< ", Total pre-calculated predecessors : " << all_init_preds
		<< ", Sent : " << (all_send_msgs/nranks)
		<< ", Received : " << (all_recv_msgs/nranks)
		<< ", Local :" << (all_local_msgs/nranks)
		<< ", Predecessors in Set Intersection : " << (all_preds/nranks)
		<< ", Successors in Set Intersection : " << (all_succs/nranks)
		<< ", Predecessors of successors in Set Intersection : " << (all_succ_preds/nranks)
		<< ", Average set intersection time : " << (tot_avg_set_int_times/total_set_inters) << " nanoseconds"
		<< ", Total set intersections : " << (total_set_inters/nranks)
		<< ", Total saved msgs : " << (all_saved_msgs/nranks)
		<< ", Total set comparisons : " << all_comparisons
		<< ", All bytes sent over : " << all_bytes_over_nw 
		<< ", Bytes sent per rank :" << (all_bytes_over_nw/nranks)
		<< ", With Ordering ? : " << with_ordering
		<< ", Predicted PSP Bytes : " << (all_predicted_psp_bytes * sizeof(Vertex))
		<< ", Predicted SPS Bytes : " << (all_predicted_sps_bytes * sizeof(Vertex))
		<< ", Max predecessor block : " << vmax_pred_block
		<< ", Max successor (group)block : " << vmax_succ_block
		<< std::endl;
    }

  }
#endif


private:
  //  template<typename Vertex>
  struct block_msg {
  private:
    typename std::vector<Vertex>::iterator sucbegin;
    uint64_t succount;
    typename std::vector<Vertex>::iterator predbegin;
    uint64_t predcount;
    Vertex* successor_array;
    Vertex* predecessor_array;

  public:
    block_msg(typename std::vector<Vertex>::iterator s,
	      uint64_t sc,
	      typename std::vector<Vertex>::iterator p,
	      uint64_t pc): sucbegin(s),
			 succount(sc),
			 predbegin(p),
			 predcount(pc),
			 predecessor_array(NULL),
			 successor_array(NULL){}

    block_msg(): succount(0), predcount(0), 
		 predecessor_array(NULL),
		 successor_array(NULL){}

    block_msg(const block_msg& bmsg): sucbegin(bmsg.sucbegin),
				      succount(bmsg.succount),
				      predbegin(bmsg.predbegin),
				      predcount(bmsg.predcount),
				      successor_array(bmsg.successor_array),
				      predecessor_array(bmsg.predecessor_array){}


    size_t get_size() {       
      // additional 1 is to store the number of successors
      return (sizeof(Vertex)*(succount+predcount+1/*number of successors*/)); 
    }

    void serialize(COALESCE_TYPE* buf) {
      Vertex* vbuf = (Vertex*)buf;
      vbuf[0] = (Vertex)succount; // Vertex is also uint64_t so should not be an issue
      // copy successors
      std::copy(sucbegin, sucbegin+succount, vbuf+1);
      // copy predecessors to the buffer
      std::copy(predbegin, predbegin+predcount, vbuf+succount+1);
    }

    void deserialize(COALESCE_TYPE* buf, uint64_t bytecount) {
      Vertex* array = (Vertex*)buf;
      succount = array[0];
      assert(succount != 0);
      successor_array = array+1;
      predcount = (bytecount/sizeof(Vertex)) - (1+succount);
      assert(predcount != 0);
      predecessor_array = array+(1+succount);
    }

    Vertex* get_predecessor_array() {
      assert(predecessor_array != NULL);
      return predecessor_array;
    }

    Vertex* get_successor_array() {
      assert(successor_array != NULL);
      return successor_array;
    }

    uint64_t get_pred_count() {
      return predcount;
    }

    uint64_t get_suc_count() {
      return succount;
    }
  };


  void initialize();

  void process(int tid, Vertex s, const Vertex* buff, const int count);


  template<typename ItePred, typename IteSuc>
  void many_set_intersection(ItePred itepbegin, 
			     uint64_t predcount,
			     IteSuc itesbegin,
			     uint64_t succount,
			     int tid,
			     bool local=false);
				 
				 
  void populate_succs_preds(const work_item_t& wi);

  template<typename SizeType>
  inline SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }

  struct vertex_logicalid_comparator {
  private:
    triangle_counting_blocked* self;
  public:
    vertex_logicalid_comparator(triangle_counting_blocked& s) : self(&s) {}

    bool operator()(Vertex v1, Vertex v2) const {
      return (self->logical_id(v1) > self->logical_id(v2));
    }
  };

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


#ifdef TC_STATS
  std::vector<uint64_t> send_msgs;
  std::vector<uint64_t> recv_msgs;
  std::vector<uint64_t> local_msgs;
  std::vector<uint64_t> total_preds;
  std::vector<uint64_t> total_succs;
  std::vector<uint64_t> total_succs_preds;
  std::vector<uint64_t> total_setints;
  std::vector<uint64_t> tot_set_int_time;
  std::vector<uint64_t> saved_msgs;
  std::vector<uint64_t> tot_comparisons;
  std::vector<uint64_t> num_bytes_sent_over_nw;
  std::vector<uint64_t> init_succs;
  std::vector<uint64_t> init_preds;
  std::vector<uint64_t> predicted_sps_bytes;
  std::vector<uint64_t> predicted_psp_bytes;
  std::vector<uint64_t> max_pred_block;
  std::vector<uint64_t> max_succ_block;
#endif


private:
  const int dummy_first_member_for_init_order;
  const Graph& g;
  amplusplus::transport& transport;
  const int nthreads;
  const IdDistribution& id_distribution;
  const OwnerMap& owner; 
  int core_offset;
  NeighborMap& vertex_predecessors;
  NeighborMap& vertex_successors;

  uint64_t block_size;
  uint64_t suc_block_size;
  uint64_t coalescing_size = 1 << 21;
  // AM++ message type
  tm_type msg_type;
  RelaxMessage relax_msg;

  shared_ptr<amplusplus::detail::barrier> t_bar;

  time_type start_time;
  time_type end_time;

  std::map<Vertex, spinlock*> suc_locks;
  std::map<Vertex, spinlock*> pred_locks;

#ifdef TRIANGLE_ENUMERATE
  work_item_t** all_triangles;
#endif
  uint64_t* threaded_triangle_indexes;

#ifdef PRINT_DEBUG
  uint64_t no_sends = 0;
  uint64_t no_receives = 0;
  Vertex lastsucc;
  uint64_t lastcount = 0;
#endif

};

#define TC_PARAMS                                   \
      typename Graph, typename IdDistribution, typename NeighborMap, typename MessageGenerator

#define TC_TYPE                                    \
      triangle_counting_blocked<Graph, IdDistribution, NeighborMap, MessageGenerator>


#ifdef TRIANGLE_ENUMERATE
#define MAX_TRIANGLE_COUNT 100000000
#endif

template<TC_PARAMS>
void
TC_TYPE::initialize() {

  //  relax_msg.set_handler(processing_function(*this));

  relax_msg.set_handler(degree_processing_function(*this));
  msg_type.set_handler(processing_function(*this));

  //  threaded_triangles.resize(nthreads);
  threaded_triangle_indexes = new uint64_t[nthreads];

  for (int i=0; i < nthreads; ++i) {
    threaded_triangle_indexes[i] = 0;
  }

  // check the coalescing size is at least upto the block size
  if (coalescing_size < block_size) {
    std::cout << "[WARNING] Coalescing size must be greater than the block size. " << std::endl;
    assert(false);
  }

  BGL_FORALL_VERTICES_T(v, g, Graph) {
    suc_locks.insert(std::make_pair(v, new spinlock()));
    pred_locks.insert(std::make_pair(v, new spinlock()));
  }

#ifdef TC_STATS
  send_msgs.resize(nthreads, 0);
  recv_msgs.resize(nthreads, 0);
  local_msgs.resize(nthreads, 0);
  total_preds.resize(nthreads, 0);
  total_succs.resize(nthreads, 0);
  total_succs_preds.resize(nthreads, 0);
  total_setints.resize(nthreads, 0);
  tot_set_int_time.resize(nthreads, 0);
  saved_msgs.resize(nthreads, 0);
  tot_comparisons.resize(nthreads, 0);
  num_bytes_sent_over_nw.resize(nthreads, 0);
  init_succs.resize(nthreads, 0);
  init_preds.resize(nthreads, 0);
  predicted_sps_bytes.resize(nthreads, 0);
  predicted_psp_bytes.resize(nthreads, 0);
  max_pred_block.resize(nthreads, 0);
  max_succ_block.resize(nthreads, 0);
#endif

#ifdef TRIANGLE_ENUMERATE
  all_triangles = new work_item_t*[nthreads];
#endif
}


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
template<typename ItePred, typename IteSuc>
void
TC_TYPE::many_set_intersection(ItePred itepbegin,
			       uint64_t predcount,
			       IteSuc itesbegin,
			       uint64_t succount,
			       int tid,
			       bool local) {


#ifdef TC_STATS
  total_preds[tid] += predcount;
  total_succs[tid] += succount;
#endif

  ItePred itepend = itepbegin + predcount;

  typedef typename std::vector<Vertex>::iterator SucVertexIter_t;
  typedef typename std::pair<SucVertexIter_t,
			     SucVertexIter_t> IteratorPair_t;

  std::vector<IteratorPair_t> sucpreds(succount);
  for (int i=0; i < succount; ++i) {
    Vertex s = *(itesbegin+i);
    sucpreds[i].first = vertex_predecessors[s].begin();
    sucpreds[i].second = vertex_predecessors[s].end();
#ifdef TC_STATS
    total_succs_preds[tid] += vertex_predecessors[s].size();
#endif
  }


  for(int i=0; i < sucpreds.size(); ++i) {

    counting_output_iterator output_ite;
    output_ite = std::set_intersection(itepbegin, itepend,
				       sucpreds[i].first,
				       sucpreds[i].second,
				       output_ite);

    threaded_triangle_indexes[tid] += output_ite.get_count();;

  }
}

template<TC_PARAMS>
void
TC_TYPE::populate_succs_preds(const work_item_t& wi) {
  Vertex d = wi.first;
  Vertex s = wi.second.first;
  uint32_t deg = wi.second.second;

  if (deg < out_degree(d, g)) {
    suc_locks[d]->lock();
    vertex_successors[d].push_back(s);
    suc_locks[d]->unlock();
  } else if (deg > out_degree(d, g)) {
    pred_locks[d]->lock();
    vertex_predecessors[d].push_back(s);
    pred_locks[d]->unlock();
  } else {
    if (logical_id(d) > logical_id(s)) {
      suc_locks[d]->lock();
      vertex_successors[d].push_back(s);
      suc_locks[d]->unlock();
    } else {
      pred_locks[d]->lock();
      vertex_predecessors[d].push_back(s);
      pred_locks[d]->unlock();
    }
  }

}


template<TC_PARAMS>
void
TC_TYPE::process(int tid, Vertex s, const Vertex* buff, const int count) {

  counting_output_iterator output_ite;

  auto begin = vertex_predecessors[s].begin();
  auto end = vertex_predecessors[s].end();

  output_ite = std::set_intersection(buff, buff+count,
				     begin,
				     end,
				     output_ite);

  (threaded_triangle_indexes[tid]) += output_ite.get_count();
}


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

#ifdef CRAYPAT
    if (PAT_region_begin ( 1, "tcrun" ) == PAT_API_FAIL) {
      std::cout << "PAT begin failed ! " << std::endl;
      assert(false);
    }
#endif

    time_type atcalps = get_time();

#ifdef TC_NO_ORDERING
    // Parallel iterate over all the vertices and collect predecessors and successors
    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	Vertex u = target(e, g);
	if (logical_id(v) > logical_id(u)) { // predecessors
	  vertex_predecessors[v].push_back(u);
	} else {
	  vertex_successors[v].push_back(u);
	}
      }

      vertex_predecessors[v].shrink_to_fit();
      vertex_successors[v].shrink_to_fit();

      std::sort(vertex_predecessors[v].begin(), vertex_predecessors[v].end());
      std::sort(vertex_successors[v].begin(), vertex_successors[v].end());

#ifdef TC_STATS
      // calculate the bytes transferred if we do psp or sps
      size_t p = vertex_predecessors[v].size();
      size_t s = vertex_successors[v].size();

      std::set<amplusplus::transport::rank_type> setranks;

      if (p > max_pred_block[tid])
	max_pred_block[tid] = p;

      // calculate the predecessor ranks and successor ranks
      int pranks = 0;
      amplusplus::transport::rank_type last_rank = -1;
      uint64_t nonlocalpreds = 0;
      for (auto i = 0; i < p; ++i) {
	auto ite = vertex_predecessors[v].begin() + i;
	amplusplus::transport::rank_type dest = get(owner, (*ite));
	if (dest != transport.rank()) {
	  setranks.insert(dest);
	  ++nonlocalpreds;
	}
      }

      pranks = setranks.size();
      setranks.clear();

      uint64_t groupsz = 0;
      int sranks = 0;
      uint64_t nonlocalsuccs = 0;
      last_rank = -1;
      for (auto i = 0; i < s; ++i) {
	auto ite = vertex_successors[v].begin() + i;
	amplusplus::transport::rank_type dest = get(owner, (*ite));
	if (dest != transport.rank()) {
	  if (last_rank != dest) {
	    last_rank = dest;
	    if (groupsz > max_succ_block[tid])
	      max_succ_block[tid] = groupsz;

	    groupsz = 0;

	  }
	  setranks.insert(dest);
	  ++groupsz;
	  ++nonlocalsuccs;
	}
      }

      sranks = setranks.size();
      setranks.clear();

      // suppose we are doing psp and calculate the number bytes transferred
      if ((p != 0) && (s != 0)) {
	auto pspbytes = (p*sranks) + nonlocalsuccs;
	auto spsbytes = (s*pranks) + nonlocalpreds;

	predicted_psp_bytes[tid] += pspbytes;
	predicted_sps_bytes[tid] += spsbytes;
      }


      init_succs[tid] += vertex_successors[v].size();
      init_preds[tid] += vertex_predecessors[v].size();
#endif

    }
#else
    { 
      amplusplus::scoped_epoch epoch(transport); 

      BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  Vertex u = target(e, g);
	  work_item_t wi = construct_wi(u, v, out_degree(v, g));
	  relax_msg.send(wi);
	}
      }
    }

    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
      std::sort(vertex_predecessors[v].begin(), vertex_predecessors[v].end());
      std::sort(vertex_successors[v].begin(), vertex_successors[v].end());

#ifdef TC_STATS
      // calculate the bytes transferred if we do psp or sps
      size_t p = vertex_predecessors[v].size();
      size_t s = vertex_successors[v].size();

      // calculate the predecessor ranks and successor ranks
      int pranks = 0;
      uint64_t nonlocalpreds = 0;
      std::set<amplusplus::transport::rank_type> setranks;

      if (p > max_pred_block[tid])
	max_pred_block[tid] = p;

      amplusplus::transport::rank_type last_rank = -1;
      for (auto i = 0; i < p; ++i) {
	auto ite = vertex_predecessors[v].begin() + i;
	amplusplus::transport::rank_type dest = get(owner, (*ite));
	if (dest != transport.rank()) {
	  setranks.insert(dest);
	  ++nonlocalpreds;
	}
      }

      pranks = setranks.size();
      setranks.clear();
      
      uint64_t groupsz = 0;
      int sranks = 0;
      uint64_t nonlocalsuccs = 0;
      for (auto i = 0; i < s; ++i) {
	auto ite = vertex_successors[v].begin() + i;
	amplusplus::transport::rank_type dest = get(owner, (*ite));
	if (dest != transport.rank()) {
	  if (last_rank != dest) {
	    last_rank = dest;
	    if (groupsz > max_succ_block[tid])
	      max_succ_block[tid] = groupsz;

	    groupsz = 0;
	  }
	  setranks.insert(dest);
	  ++groupsz;
	  ++nonlocalsuccs;
	}
      }

      sranks = setranks.size();
      setranks.clear();

      // suppose we are doing psp and calculate the number bytes transferred
      if ((p != 0) && (s != 0)) {
	auto pspbytes = (p*sranks) + nonlocalsuccs;
	auto spsbytes = (s*pranks) + nonlocalpreds;

	predicted_psp_bytes[tid] += pspbytes;
	predicted_sps_bytes[tid] += spsbytes;
      }

      init_succs[tid] += vertex_successors[v].size();
      init_preds[tid] += vertex_predecessors[v].size();
#endif

    }
#endif

    // debugging -- remove later
    /*if (tid == 0)
      write_to_file();*/

    t_bar->wait();
    
    // debugging -- remove
    //{ amplusplus::scoped_epoch epoch(transport); }

#ifdef PRINT_DEBUG
    t_bar->wait();
    if (tid == 0) {
      size_t cap = 0;
      size_t totsz = 0;
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	cap += vertex_predecessors[v].capacity();
	cap += vertex_successors[v].capacity();
	totsz += vertex_predecessors[v].size();
	totsz += vertex_successors[v].size();
      }

      std::cout << "Rank : " << transport.rank() << " capacity : " << cap << " size : " << totsz << std::endl;
    }        
#endif
    // For debugging
    { amplusplus::scoped_epoch epoch(transport); }

    typedef std::pair<typename std::vector<Vertex>::iterator, uint64_t> SuccIterSizePair_t;
    std::vector< SuccIterSizePair_t > rank_successors(transport.size());

    for (amplusplus::transport::rank_type r=0; 
	 r < transport.size(); ++r) {
      rank_successors[r].second = 0;
    }

    time_type etcalps = get_time();
    if ((transport.rank()==0) && (tid == 0))
      std::cout << "Time to calculate predecessors and successors : " << (etcalps-atcalps) << std::endl;

    uint64_t offset = 0;
    // should come before begin epoch
    start_time = get_time();
    // Start the algorithm

    t_bar->wait();

    std::pair<typename boost::graph_traits<Graph>::vertex_iterator,
    	      typename boost::graph_traits<Graph>::vertex_iterator> itepair = vertices(g);

    typename boost::graph_traits<Graph>::vertex_iterator startite = itepair.first;

    time_type atint = get_time();
    {
      amplusplus::scoped_epoch epoch(transport);
      for(; startite != itepair.second; ++startite) {

	// no successors or predecessors, then continue
	if ((vertex_successors[*startite].size() == 0) ||
	    (vertex_predecessors[*startite].size() == 0))
	  continue;

	offset = (local_id(*startite) + tid) % nthreads;

#ifdef PRINT_DEBUG
	if (transport.rank() == 0) {
	  if (local_id(*startite) % 1000) {
	    std::cout << "R: " << transport.rank() << "vertex : " << *startite << std::endl;
	    std::cout << "R: " << transport.rank() << "offset : " << offset << std::endl;	
	  }
	}
#endif

	uint64_t pred_count = vertex_predecessors[*startite].size();
	uint64_t succ_count = vertex_successors[*startite].size();

#ifdef PRINT_DEBUG
	std::cout << "pred_count : " << pred_count << std::endl;
	std::cout << "succ_count : " << succ_count << std::endl;
#endif
	//	for (auto pos = offset; pos < pred_count; pos=pos+nthreads) {
	auto numblocks = (pred_count+block_size-1)/block_size;
	auto numsucblocks = (succ_count+suc_block_size-1)/suc_block_size;
	
	//std::cout << "pred blocks : " << numblocks << std::endl;
	//std::cout << "suc blocks : " << numsucblocks << std::endl;
	
	for (auto pos = offset; pos < numblocks; pos=pos+nthreads) {
	  auto begin = pos * block_size;
	  size_t size = std::min(block_size, (pred_count-(pos*block_size)));
	  auto end = begin + size;

	  auto predbegin = vertex_predecessors[*startite].begin() + begin;
	  auto predend = vertex_predecessors[*startite].begin() + end;

	  for (uint64_t sucblkpos = 0; sucblkpos < numsucblocks; ++sucblkpos) {
#ifdef PRINT_DEBUG
	    std::cout << "begin :" << begin << " end : " << end << std::endl;
#endif
	    auto sucpos = sucblkpos * suc_block_size;
	    auto sucposend = sucpos + std::min(suc_block_size, (succ_count-(sucblkpos*suc_block_size)));
	    
	    for(; sucpos < sucposend; ++sucpos) {

	      Vertex send_succ = *(vertex_successors[*startite].begin() + sucpos);
	      amplusplus::transport::rank_type dest = get(owner, send_succ);	    
	      if (rank_successors[dest].second == 0) {
		rank_successors[dest].first = (vertex_successors[*startite].begin() + sucpos);
	      }
	      
	      rank_successors[dest].second++;
	    }

	    // now send these to appropriate ranks
	    // first send to remote ranks
	    for (amplusplus::transport::rank_type r=0; 
		 r < transport.size(); ++r) {
	      if ((r != transport.rank()) && (rank_successors[r].second != 0)) {

#ifdef TC_STATS
		send_msgs[tid]++;
#endif

		block_msg bmsg(rank_successors[r].first, rank_successors[r].second,
			       predbegin, size);
#ifdef TC_STATS
		num_bytes_sent_over_nw[tid] += (bmsg.get_size() - sizeof(Vertex));
#endif

		msg_type.send(bmsg, r);
		rank_successors[r].second = 0;
	      }
	    }

	    // now do the local rank
	    if (rank_successors[transport.rank()].second != 0) {
#ifdef TC_STATS
	      local_msgs[tid]++;
	      auto sistart = std::chrono::system_clock::now();
#endif


	      many_set_intersection(predbegin, size,
				    rank_successors[transport.rank()].first,
				    rank_successors[transport.rank()].second,
				    tid, true /*local call*/);
#ifdef TC_STATS
	      auto siend = std::chrono::system_clock::now();
	      total_setints[tid]++;
	      auto durationforthis 
		= std::chrono::duration_cast<std::chrono::nanoseconds>(siend - sistart).count();
	      tot_set_int_time[tid] += durationforthis;
#endif

	      rank_successors[transport.rank()].second = 0;
	    }
	  }
	}
      }
    }

    time_type etint = get_time();
    if (tid == 0)
      std::cout << "Set intersection time : " << (etint-atint) << std::endl;
  

#ifdef PRINT_DEBUG
    std::cout << "End of epoch ............................." << std::endl;
    std::cout << "Sends : " << no_sends 
	      << ", Receives : " << no_receives 
	      << ", Last Succ : " << lastsucc 
	      << ", Last Count : " << lastcount << std::endl;

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
  processing_function(triangle_counting_blocked& self) : self(&self) {}
  
  void operator() (const amplusplus::transport::rank_type src, 
		   block_msg& msg) const {
    int tid = amplusplus::detail::get_thread_id();

#ifdef TC_STATS
    self->recv_msgs[tid]++;
#endif

#ifdef TC_STATS
    auto sistart = std::chrono::system_clock::now();
#endif

    self->many_set_intersection(msg.get_predecessor_array(),
				msg.get_pred_count(),
				msg.get_successor_array(),
				msg.get_suc_count(),
    				tid);

#ifdef TC_STATS
    auto siend = std::chrono::system_clock::now();
    self->total_setints[tid]++;
    auto currenttot = self->tot_set_int_time[tid] * (self->total_setints[tid]-1) +
      std::chrono::duration_cast<std::chrono::milliseconds>(siend - sistart).count();
    self->tot_set_int_time[tid] = currenttot / self->total_setints[tid];
#endif

  }

protected:
  triangle_counting_blocked* self;
};

template<TC_PARAMS>
struct TC_TYPE::
degree_processing_function {
  
  degree_processing_function() : self(NULL) {}
  degree_processing_function(triangle_counting_blocked& self) : self(&self) {}
  
  void operator() (const work_item_t& data) const {
    int tid = amplusplus::detail::get_thread_id();
    self->populate_succs_preds(data);
  }

protected:
  triangle_counting_blocked* self;
};


}}}
#endif
