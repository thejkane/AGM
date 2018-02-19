#ifndef __AMPP_RUNTIME_HPP__
#define __AMPP_RUNTIME_HPP__

#include "runtime.hpp"
#include <boost/graph/agm/runtime/runtime_common.hpp>
#include <boost/graph/agm/util/runtime_stats.hpp>
#include <boost/graph/util/utils.hpp>

#include <am++/make_mpi_datatype.hpp>
#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>
#include <unistd.h>

// NUMA STUFF
#ifndef __APPLE__
#include <numaif.h>
#include <numa.h>
#endif

namespace boost { namespace graph { namespace agm {

// TODO encapsulate stat function calls inside PBGL2_WORK_STATS macro
template<typename work_item,
         typename receive_handler,         
         typename post_order_processing_function,
         typename owner_map,
         typename message_generator,
         typename ampp_runtime_gen>
class ampp_runtime;

template<typename MessageGenerator,
         typename OwnerMap>
struct ampp_runtime_gen {
public:
  amplusplus::transport& transport;
  MessageGenerator& msg_gen;
  OwnerMap& owner;
  runtime_stats& rt_stats;  
  int core_offset;
  instance_params& rtparams;
  
  ampp_runtime_gen(amplusplus::transport& _t,
                   MessageGenerator& _mgen,
                   OwnerMap& _owner,
                   runtime_stats& _stat,
                   int offs,
                   instance_params& _rtparams): transport(_t),
                                                msg_gen(_mgen),
                                                owner(_owner),
                                                rt_stats(_stat),
                                                core_offset(offs),
                                                rtparams(_rtparams){}
                                             
  
  template <typename WorkItem,
            typename ReceiveHandler,
            typename PostOrderProcessingFunction>
  struct inner {
    typedef ampp_runtime<WorkItem,
                         ReceiveHandler,
                         PostOrderProcessingFunction,
                         OwnerMap,
                         MessageGenerator,
                         ampp_runtime_gen<MessageGenerator,
                                          OwnerMap> > type;
  };
};

template<typename WorkItem,
         typename RelaxMessage>
struct out_sender {
private:
  RelaxMessage& relax_msg;

public:
  out_sender(RelaxMessage& _relax_msg):relax_msg(_relax_msg){}
  inline void push(const WorkItem& wi, int tid) {
    relax_msg.send(wi);
  }
};
      
template<typename work_item,
         typename receive_handler,         
         typename post_order_processing_function,
         typename owner_map,
         typename message_generator,
         typename ampp_runtime_gen>
class ampp_runtime : public runtime_base {

public:  
  typedef work_item work_item_t;
  
  typedef typename message_generator::template call_result<work_item, 
                                                           receive_handler,
                                                           work_item_owner<owner_map, work_item>,
                                                           amplusplus::no_reduction_t>::type
  RelaxMessage;

  typedef out_sender<work_item, RelaxMessage> OutSenderType;
  
private:
  int dummy;
  amplusplus::transport& transport;
  MPI_Datatype dt;
  owner_map& owner;  
  int core_offset;
  RelaxMessage relax_msg;  
  runtime_stats& stats;
  post_order_processing_function& post_pf;
  OutSenderType out_set;
  boost::shared_ptr<amplusplus::detail::barrier> t_bar;
  std::vector<int> thread_numa_map; // maps thread id to numa domain
  std::vector<int> numa_thread_indexes; // finds the thread index within a numa domain
  std::vector<int> numa_main_thread_map; // maps numa domain to its main thread
  std::vector<int> numa_thread_count_map;
  std::vector< boost::shared_ptr<amplusplus::detail::barrier> > numa_thread_barriers;
  

  void initialize() {
    // initialize the thread barrier
    t_bar.reset(new amplusplus::detail::barrier(this->threads));
    numa_init();
  }

  void numa_init() {
#ifndef __APPLE__
    // check whether NUMA available
    if (numa_available() < 0) {
      warn("[WARN]System does not support numa.");
      this->numa_domains = 1;
      return;
    }

    // NUMA is available
    this->numa_domains = (numa_max_node() + 1);
    info("Number of NUMA domains : ", this->numa_domains);
#else
    this->numa_domains = 1;
#endif

    thread_numa_map.resize(this->threads, 0);    
    numa_thread_indexes.resize(this->threads);
    numa_main_thread_map.resize(this->numa_domains, std::numeric_limits<int>::max());
    numa_thread_count_map.resize(this->numa_domains, 0);

    numa_thread_barriers.resize(this->numa_domains);
  }
  
public:
  ampp_runtime(const ampp_runtime_gen& _rgen,
               post_order_processing_function& _sendpf):
    runtime_base(_rgen.rtparams.threads,
                 1,
                 _rgen.transport.size()),
    dummy((amplusplus::register_mpi_datatype<work_item>(), 0)),
    transport(_rgen.transport),
    dt(amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item>())),
    owner(_rgen.owner),
    core_offset(_rgen.core_offset),
    relax_msg(_rgen.msg_gen,
              transport,
              work_item_owner<owner_map, work_item>(owner),
              amplusplus::no_reduction),    
    stats(_rgen.rt_stats),
    post_pf(_sendpf),
    out_set(relax_msg){

    initialize();
  }

  void register_receiver(receive_handler& _handler) {
    relax_msg.set_handler(_handler);
  }

  void print_numa_data() {
    info("Number of threads : ", this->threads);
    info("Number of NUMA domains : ", this->numa_domains);
    for (int i=0; i < numa_thread_count_map.size(); ++i) {
      info("NUMA id : " , i, ", threads for NUMA node : ", numa_thread_count_map[i]);
    }

    for (int i=0; i < numa_main_thread_map.size(); ++i) {
      info("NUMA id : " , i, ", main thread id : ", numa_main_thread_map[i]);
    }

    for(int i=0; i < thread_numa_map.size(); ++i) {
      info("Thread id : ", i, ", NUMA node : ", thread_numa_map[i]);
    }

    for(int i=0; i < numa_thread_indexes.size(); ++i) {
      info("Thread id : ", i, ", Index within NUMA : ", numa_thread_indexes[i]);
    }
  }

  void initialize_per_thread(int tid) {
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

#ifndef __APPLE__    
    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    thread_numa_map[tid] = node;
    ++numa_thread_count_map[node];
    if (tid < numa_main_thread_map[node])
      numa_main_thread_map[node] = tid;
#else
    numa_main_thread_map[0] = 0;
    ++numa_thread_count_map[0];
#endif    

    // this is a map from tid to to index within the numa domain
    numa_thread_indexes[tid] = tid;

    t_bar->wait();
    if (tid == 0) {
      for(int i=0; i < this->numa_domains; ++i) {
	int numaindex = 0;
	for(int j=0; j < numa_thread_indexes.size(); ++j) {
	  if (thread_numa_map[j] == i){
	    numa_thread_indexes[j] = numaindex;
	    ++numaindex;
	  }
	}

	info("Initializing barrier for NUMA domain :", i, ", to : ", numaindex);
	// if no threads specified to run in numa domain numaindex, leave
	// the barrier null
	if (numaindex != 0) {
	  numa_thread_barriers[i].reset(new amplusplus::detail::barrier(numaindex));
	}
      }

      print_numa_data();
    }
    t_bar->wait();
  }

  void set_transport_threads(int t) {
    transport.set_nthreads(t);
  }

  void send(const work_item& wi, int tid) {
    post_pf(wi, tid, out_set);
    stats.increment_sends(tid);    
    this->decrease_activity_count(tid, 1);
  }

  inline int find_numa_node(int tid) {
    return thread_numa_map[tid];
  }
  
  void do_all_gather(void* _ptosend,
                     int sendcnt,
                     void* _precv,
                     int recvcnt) {
    
    stats.increment_all_to_alls();
    MPI_Allgather(_ptosend, sendcnt, dt, _precv, recvcnt, dt, MPI_COMM_WORLD);
  }

  int get_nthreads_for_numa_domain(int tid) {
    int numa = find_numa_node(tid);    
    return numa_thread_count_map[numa]; // TODO change
  }

  int get_thread_index_in_numa_domain(int tid) {
    return numa_thread_indexes[tid]; 
  }

  void wait_for_numa_domain_threads_to_reach_here(int tid) {
    int numa = find_numa_node(tid);
    numa_thread_barriers[numa]->wait();
  }

  bool is_main_thread_in_numa_domain(int tid) {
    int numa = find_numa_node(tid);
    return (numa_main_thread_map[numa]==tid);
  }
  
  bool begin_epoch(int tid){
    stats.increment_epochs(tid);
    transport.begin_epoch();
  }

  void end_epoch() {
    transport.end_epoch();
  }

  unsigned long end_epoch_with_value(const unsigned long& _read_value) {
    return transport.end_epoch_with_value(_read_value);
  }

  void pull_work(int tid) {
    // TODO
  }
  
  void synchronize(){
    amplusplus::scoped_epoch epoch(transport); 
  }

  void wait_for_threads_to_reach_here(int tid){
    t_bar->wait();
  }

  bool is_main_thread(int tid) {
    return (tid == 0);
  }    
};

}}}
#endif
