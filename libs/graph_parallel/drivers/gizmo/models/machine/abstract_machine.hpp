#ifndef __AGM_DISTRIBUTED_ABSTRACT__
#define __AGM_DISTRIBUTED_ABSTRACT__

#include <cstdint>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

/*class abstract_machine_base {
public:
  virtual void flush(uint64_t t) = 0;
  };*/

//TODO -- remove this
void invoke_gizmo_cb(void* _m, uint64_t t) {
  //assert(_m != NULL);
  //  abstract_machine_base* pBase = (abstract_machine_base*)_m;
  //pBase->flush(t);
}

template<typename WorkItem,
	 typename StrictWeakOrderingRelation,
	 typename ProcessingFunction,
	 typename GraphDistribution>
class abstract_machine {

public:
  struct Comp {
  public:
    bool operator()(const WorkItem& left, const WorkItem& right) const {
      return (std::get<2>(left) > std::get<2>(right));
    }
  };

  //  typedef std::priority_queue<WorkItem, std::vector<WorkItem>, Comp> PriorityQueueType;
  typedef std::queue<WorkItem> PriorityQueueType;
  typedef buckets<WorkItem, StrictWeakOrderingRelation> buckets_t;


  class simulating_worker {
  private:
    uint64_t ts;
    int tid;
    std::queue<WorkItem> receive_q;
    int flow;
    ProcessingFunction& pf;
    buckets_t* pall_buckets;
    static int flow_control;

  public:
    simulating_worker(int _tid, ProcessingFunction& _pf,
		      buckets_t* _pall): ts(0), 
					 tid(_tid),
					 flow(0), 
					 pf(_pf),
					 pall_buckets(_pall) {}

    static void set_flow_control(int f) {
      flow_control = f;
    }

    void flush() {
      // flush the receive queue
      while(!receive_q.empty()) {
	WorkItem wi = receive_q.front();
	receive_q.pop();
	if (std::get<2>(wi) > ts) {
	  ts = std::get<2>(wi);
	}

	std::vector<WorkItem> out;
	pf(wi, out, ts);
	if (!out.empty()) {
	  push_work_items(out.begin(), out.end());
	}
      }
    }

    int get_tid() {
      return tid;
    }

    void print_ts() {
      std::cout << "TID : " << tid << " Timestamp : " << ts << std::endl;
    }


    template<typename ContainerIterator>
    void push_work_items(ContainerIterator begin, ContainerIterator end) {
      pall_buckets->push(begin, end);
    }

    uint64_t get_timestamp() const {
      return ts;
    }

    void push(WorkItem& wi) {
      receive_q.push(wi);
    }

    bool empty() {
      return receive_q.empty();
    }

    template<typename queue_t>
    void send(WorkItem& wi, 
	      uint64_t totlatency /*pop_overhead + send_overhead + transmission_time*/,
	      queue_t& outgoing) {
      //assert(std::get<2>(wi) <= ts);
      ts += totlatency;
      std::get<2>(wi) = ts;
      //    debug("Pushing to work items to rank queueus ...");
      ++flow;
      outgoing.push(wi);
    
      if (flow >= flow_control) {
	debug("Invoking flow control ...");
	flow = 0;
	flush();
      }      
    }
  };

  struct rank_data {
  private:
    int num_workers;
    std::vector<simulating_worker*> workers;
    boost::random::mt19937 rng;

    /*    struct TimeStampComp {
    public:
      bool operator()(const simulating_worker& w1,
		      const simulating_worker& w2) const {
	return (w1.get_timestamp() < w2.get_timestamp());
      }
      };*/

  public:
    rank_data(int _w, ProcessingFunction& _pf,
	      buckets_t* _pall) : num_workers(_w) {
      for(int i=0; i < num_workers; ++i) {
	workers.push_back(new simulating_worker(i, _pf, _pall));
      }
    } 

    ~rank_data() {
      for(int i=0; i < num_workers; ++i) {
	delete workers[i];
      }
      workers.clear();
    }

    void print_ts() {
      for (int i=0; i < num_workers; ++i) {
	workers[i]->print_ts();
      }
    }

    uint64_t get_highest_ts() {
      uint64_t highest = 0;
      for (int i=0; i < num_workers; ++i) {
	if (workers[i]->get_timestamp() > highest) {
	  highest = workers[i]->get_timestamp();
	}
      }

      return highest;
    }

    uint64_t get_lowest_ts() {
      uint64_t lowest = UINT64_MAX;
      for (int i=0; i < num_workers; ++i) {
	if (workers[i]->get_timestamp() < lowest) {
	  lowest = workers[i]->get_timestamp();
	}
      }

      return lowest;
    }

    static bool time_stamp_comp(const simulating_worker* pw1,
			 const simulating_worker* pw2) {
      return (pw1->get_timestamp() < pw2->get_timestamp());
    }

    // TODO generalize this
    simulating_worker& get_smallest_ts_worker(uint64_t ts) {
      assert(workers.size() > 0);
      // do a sort and get the smallest ts worker
      //      std::sort(workers.begin(), workers.end(), TimeStampComp());
      std::sort(workers.begin(), workers.end(), time_stamp_comp);
      return (*workers[0]);
    }

    simulating_worker& get_random_worker() {
      boost::random::uniform_int_distribution<> rnddist(0,(num_workers-1));
      int i = rnddist(rng);
      assert((0 <= i) && (i < num_workers));
      return (*workers[i]);
    }

    simulating_worker& get_worker(const WorkItem& wi) {
      simulating_worker& sw = get_random_worker();
      //      std::cout << "Selected : " << sw.get_tid() << std::endl;
      return sw;
    }

    void flush_all_workers() {
      for (int i=0; i < num_workers; ++i) {
	workers[i]->flush();
      }
    }

    bool empty() {
      for (int i=0; i < num_workers; ++i) {
	if (!workers[i]->empty())
	  return false;
      }
      
      return true;
    }
  };
  
private:
  //  std::vector< PriorityQueueType > rank_queues;
  std::vector< rank_data > all_rank_data;
  int ranks;
  uint64_t num_vertices;
  ProcessingFunction& pf;
  GraphDistribution dist;
  buckets_t all_buckets;
  const int flow_control = 2;
  int threads;

public:
  abstract_machine(uint64_t n,
		   ProcessingFunction& _pf) : num_vertices(n),
					      pf(_pf){}

  bool initialize (const gizmo_config& configs) {
    info("Initializing abstract machine class ...");
    ranks = configs.get_int(KEY_RANKS);
    assert(ranks != 0);

    info("Initializing threads ...");
    threads = configs.get_int(KEY_THREADS);

    info("Initializing distributions ...");
    dist.initialize(num_vertices, ranks);
    info("Base class initializing done.");

    info("Setting flow control on workers ...");
    simulating_worker::set_flow_control(flow_control);

    info("Initializing flow control queues.");
    for(int i=0; i < ranks; ++i) {
      all_rank_data.emplace_back(threads, pf,
				 &all_buckets);
    }
  }

  uint64_t execution_time() {
    uint64_t executiont = 0;
    for (int i=0; i < ranks; ++i) {
      if (all_rank_data[i].get_highest_ts() > executiont) {
	executiont = all_rank_data[i].get_highest_ts();
      }
    }

    return executiont;
  }

  uint64_t fastest_execution_time() {
    uint64_t lowest = UINT64_MAX;
    for (int i=0; i < ranks; ++i) {
      if (all_rank_data[i].get_lowest_ts() < lowest) {
	lowest = all_rank_data[i].get_lowest_ts();
      }
    }

    return lowest;
  }


  int get_ranks() {
    return ranks;
  }

  /*  int get_total_parallel_workers() {
    return (ranks * threads);
    }*/

  virtual int get_ordering_overhead() {
    return 0;
  }

  virtual int get_distribution_overhead() {
    return 0;
  }

  virtual int get_pop_overhead() {
    return 0;
  }

  virtual int get_send_overhead() = 0;
  virtual int get_transmission_time(int source, int destination) = 0;

  /*  void send(WorkItem& wi, int tid) {
    auto target = std::get<0>(wi);
    auto src = std::get<1>(wi);
    int destination = dist.owner(target);
    int source = dist.owner(src);

    util::increment_time(wi, 
			 (get_pop_overhead() + get_send_overhead()));

    auto total = get_transmission_time(source, destination);

    std::get<2>(wi) = std::get<2>(wi) + total; 
    //    debug("Pushing to work items to rank queueus ...");
    rank_queues[destination].push(wi);
    int flow = (flows[destination].get_flow_control(tid) + 1);
    flows[destination].set_flow_control(tid, flow);
    
    if (flow >= flow_control) {
      debug("Invoking flow control ...");
      flows[destination].set_flow_control(tid, 0);
      flush_all_queues();
    }      
    }*/

  void send(WorkItem& wi) {
    auto target = std::get<0>(wi);
    auto src = std::get<1>(wi);
    int destination = dist.owner(target);
    int source = dist.owner(src);
    uint64_t ts = std::get<2>(wi);

    rank_data& dest_rank = all_rank_data[destination];
    rank_data& src_rank = all_rank_data[source];

    simulating_worker& dest_worker = dest_rank.get_worker(wi);
    simulating_worker& src_worker = src_rank.get_worker(wi);
    
    auto totaltransmission = get_pop_overhead()
      + get_send_overhead() + get_transmission_time(source, destination);

    src_worker.send(wi, totaltransmission, dest_worker);
  }


  /*  void flush(uint64_t time) {
    debug("Flush getting called ...");
    // go through every work item and flush work items that has a timestamp
    // greater than time
    for (int i=0; i < ranks; ++i) {
      while(!rank_queues[i].empty()) {
	WorkItem wi = rank_queues[i].front();	
	//	std::cout << "Workitem timestamp : " << std::get<2>(wi) << ", current time : " << time << std::endl; 
	if (std::get<2>(wi) <= time) {
	  info("Condition satisfied ...");
	  rank_queues[i].pop();

	  std::vector<WorkItem> out;
	  pf(wi, out);
	  util::invoke_call_back();
	  if (!out.empty()) {
	    debug("Pushing workitems after a flush ...");
	    push_work_items(out.begin(), out.end());
	  }
	} else {
	  break;
	}
      }
    }
    }*/

  void flush_all_queues() {
    for (int i=0; i < ranks; ++i) {
      all_rank_data[i].flush_all_workers();
    }
  }

  bool are_all_queues_empty() {
    for (int i=0; i < ranks; ++i) {
      if (!all_rank_data[i].empty())
	return false;
    }

    return true;
  }

  template<typename ContainerIterator>
  void push_work_items(ContainerIterator begin, ContainerIterator end) {
    all_buckets.push(begin, end);
  }

  inline bool all_buckets_empty() {
    return all_buckets.empty();
  }

  inline typename buckets_t::bucket_t* get_top_bucket() {
    return all_buckets.get_top_bucket();
  }

  inline void pop_bucket() {
    assert(are_all_queues_empty() && "Rank queues are not empty but pop bucket is called.");
    all_buckets.pop_bucket();
  }

  inline virtual void synchronize() {
    // TODO
  }

};

template<typename WorkItem,
	 typename StrictWeakOrderingRelation,
	 typename ProcessingFunction,
	 typename GraphDistribution>
int abstract_machine<WorkItem, StrictWeakOrderingRelation, ProcessingFunction, GraphDistribution>::simulating_worker::flow_control = 10;

#endif
