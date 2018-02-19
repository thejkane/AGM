#include <boost/tuple/tuple.hpp>
#include <limits.h>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/model/eagm_buckets.hpp>
#include <boost/graph/agm/util/eagm_config.hpp>
#include <boost/thread/locks.hpp>
#include <am++/make_mpi_datatype.hpp>
#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/util/utils.hpp>


template<typename work_item>
class test_runtime {
private:
  const int dummy;
  int threads;  
  int numa_domains;
  std::atomic<int> receives;
  std::atomic<std::uint64_t> active_count;
  std::atomic<int> activity_cnt;
  amplusplus::transport& trans;  
  MPI_Datatype dt;
  int num_syncs;
  boost::shared_ptr<amplusplus::detail::barrier> t_bar;
  std::vector< boost::shared_ptr<amplusplus::detail::barrier> > numa_thread_barriers;
  int threads_per_numa;
  
public:

  test_runtime(int _nthreads,
               int _nnuma,
               amplusplus::transport& _trans): dummy((amplusplus::register_mpi_datatype<work_item>(), 0)),
                                               threads(_nthreads),
                                               numa_domains(_nnuma),
                                               receives(0),
                                               active_count(0),
                                               trans(_trans),
                                               dt(amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item>())),
                                               num_syncs(0){

    
    t_bar.reset(new amplusplus::detail::barrier(threads));
    numa_thread_barriers.resize(numa_domains);

    assert((threads % numa_domains) == 0); // for simplicity in tests
    
    threads_per_numa = threads / numa_domains;
    for(int i=0; i < numa_domains; ++i) {
      numa_thread_barriers[i].reset(new amplusplus::detail::barrier(threads_per_numa));
    }
  }
                            
  void send(const work_item& wi, int tid) {
    receives++;    
    std::cout << "R: (" << std::get<0>(wi) << ", "
              << std::get<1>(wi) << ")" << std::endl;
    decrease_activity_count(tid, 1);
  }

  void initialize_per_thread(int tid){}
  
  void register_receiver(std::function<void(const work_item&, int)> funcrev) {
    // set message handler
  }

  int find_numa_node(int tid) {
    return (tid / threads_per_numa);
  }

  int get_nnuma_nodes() {
    return numa_domains;
  }
  
  inline int get_nthreads() {
    return threads;
  }  

  inline int get_nranks() {
    return 1;
  }

  void do_all_gather(void* _ptosend,
                     int sendcnt,
                     void* _precv,
                     int recvcnt) {
    ++num_syncs;
    MPI_Allgather(_ptosend, sendcnt, dt, _precv, recvcnt, dt, MPI_COMM_WORLD);
  }

  void increase_activity_count(int tid, uint64_t v) {
    active_count.fetch_add(v);
  }
  
  void decrease_activity_count(int tid, uint64_t v) {
    active_count.fetch_sub(v);
  }

  uint64_t get_activity_count() {
    return active_count.load();
  }
  
  int get_nthreads_for_numa_domain(int tid) {
    return threads_per_numa;
  }

  int get_thread_index_in_numa_domain(int tid) {
    return (tid % threads_per_numa);
  }

  void wait_for_numa_domain_threads_to_reach_here(int tid) {
    // numa domain
    int d = tid / threads_per_numa;
    numa_thread_barriers[d]->wait();
  }

  bool is_main_thread_in_numa_domain(int tid) {
    return ((tid%threads_per_numa) == 0);
  }
  
  bool begin_epoch(int tid){
    trans.begin_epoch();
  }
  
  bool is_epoch_completed(int tid) {
    assert(false);
  }

  void end_epoch() {
    trans.end_epoch();
  }

  unsigned long end_epoch_with_value(const unsigned long& _read_value) {
    return trans.end_epoch_with_value(_read_value);
  }

  void pull_work(int tid) {}
  
  //void synchronize(){}

  void wait_for_threads_to_reach_here(int tid){
    t_bar->wait();
  }

  bool is_main_thread(int tid) {
    return (tid == 0);
  }

  //void synchronize_numa_domain_threads(int tid){}

  int get_num_receives() {
    return receives.load();
  }

  int get_num_syncs() {
    return num_syncs;
  }
};


typedef std::tuple<uint64_t, uint64_t> work_item_t;
void init_data(std::vector<work_item_t>& wis) {
  work_item_t w1(0,5);
  work_item_t w2(0,1);
  work_item_t w3(0,9);
  work_item_t w4(0,25);
  work_item_t w5(0,15);
  work_item_t w6(0,2);
  work_item_t w7(0,34);
  work_item_t w8(0,16);
  work_item_t w9(0,8);
  work_item_t w10(0,9);
  work_item_t w11(8, 45);

  wis.push_back(w1);
  wis.push_back(w2);
  wis.push_back(w3);
  wis.push_back(w4);
  wis.push_back(w5);
  wis.push_back(w6);
  wis.push_back(w7);
  wis.push_back(w8);
  wis.push_back(w9);
  wis.push_back(w10);
  wis.push_back(w11);

}


template<typename buckets>
class thread_executor {

private:
  int nthreads;
  buckets& all_bkts;
  
public:
  thread_executor(int _n,
                  buckets& _all) : nthreads(_n),
                                   all_bkts(_all){}
  
  void operator()(int numa,
                  int tid,
                  std::vector<work_item_t>& wis) {
    for (typename std::vector<work_item_t>::size_type i = tid ;
         i < wis.size(); i+= nthreads) {
      work_item_t& wi = wis[i];
      all_bkts.push(wi, tid);      
    }  
  }  
};

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
boost::graph::agm::eagm_configs<T1, T2, T3, T4> create_eagm_configs(T1 _t1,
                                                                   T2 _t2,
                                                                   T3 _t3,
                                                                   T4 _t4) {
  typedef boost::graph::agm::eagm_configs<T1, T2, T3, T4> eagm_config_t;
  eagm_config_t config(_t1,
                       _t2,
                       _t3,
                       _t4);

  return config;
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
boost::graph::agm::eagm_configs<T1, T2, T3, T4, boost::graph::agm::pq_container>
create_eagm_config_w_node_pq(T1 _t1,
                             T2 _t2,
                             T3 _t3,
                             T4 _t4) {
  typedef boost::graph::agm::eagm_configs<T1, T2, T3, T4,
                                          boost::graph::agm::pq_container> eagm_config_t;
  eagm_config_t config(_t1,
                       _t2,
                       _t3,
                       _t4);

  return config;
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
boost::graph::agm::eagm_configs<T1, T2, T3, T4,boost::graph::agm::buffer_container,
                                boost::graph::agm::pq_container>
create_eagm_config_w_numa_pq(T1 _t1,
                             T2 _t2,
                             T3 _t3,
                             T4 _t4) {
  typedef boost::graph::agm::eagm_configs<T1, T2, T3, T4,
                                          boost::graph::agm::buffer_container,
                                          boost::graph::agm::pq_container> eagm_config_t;
  eagm_config_t config(_t1,
                       _t2,
                       _t3,
                       _t4);

  return config;
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm(T1 _t1,
               T2 _t2,
               T3 _t3,
               T4 _t4,
               amplusplus::transport& trans,
               int _numsyncs=1,
               bool process=false) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_configs<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  debug("Setting threads to 1" );
  trans.set_nthreads(1);
  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(1, 1, trans);
  
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);
  
  std::vector<work_item_t> wis;
  init_data(wis);

  int tid = 0;
  amplusplus::detail::push_thread_id_obj tobj(tid);
  
  rt.begin_epoch(0);
  
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], tid);
  }

  e.print();

  std::cout << "[===============================================================]" << std::endl;
  
  if (process) {
    e.process_global_buckets(tid);
    std::cout << "Num receives : " << rt.get_num_receives() << ", wis size : " << wis.size()
              << ", num buckets (syncs) : " << rt.get_num_syncs()
              << std::endl;

    assert(_numsyncs == rt.get_num_syncs());
    assert((rt.get_num_receives() == wis.size()) && "Number of receives is not equal to sends");
  }
}


template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm_with_node_pq(T1 _t1,
                            T2 _t2,
                            T3 _t3,
                            T4 _t4,
                            amplusplus::transport& trans,
                            int _numsyncs=1,
                            bool process=false) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4, boost::graph::agm::pq_container> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_config_w_node_pq<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);
  debug("Setting threads to 1" );
  trans.set_nthreads(1);
  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(1, 1, trans);
  
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);
  
  std::vector<work_item_t> wis;
  init_data(wis);

  int tid = 0;
  amplusplus::detail::push_thread_id_obj tobj(tid);
  
  rt.begin_epoch(0);
  
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], tid);
  }

  e.print();

  std::cout << "[===============================================================]" << std::endl;
  
  if (process) {
    e.process_global_buckets(tid);
    std::cout << "Num receives : " << rt.get_num_receives() << ", wis size : " << wis.size()
              << ", num buckets (syncs) : " << rt.get_num_syncs()
              << std::endl;

    assert(_numsyncs == rt.get_num_syncs());
    assert((rt.get_num_receives() == wis.size()) && "Number of receives is not equal to sends");
  }
}


template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm_with_numa_pq(T1 _t1,
                            T2 _t2,
                            T3 _t3,
                            T4 _t4,
                            amplusplus::transport& trans,
                            int _numsyncs=1,
                            bool process=false) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4,
                                          boost::graph::agm::buffer_container,
                                          boost::graph::agm::pq_container> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_config_w_numa_pq<T1,
                                                         T2,
                                                         T3,
                                                         T4>(_t1,
                                                             _t2,
                                                             _t3,
                                                             _t4);

  debug("Setting threads to 1" );
  trans.set_nthreads(1);
  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(1, 1, trans);
  
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);
  
  std::vector<work_item_t> wis;
  init_data(wis);

  int tid = 0;
  amplusplus::detail::push_thread_id_obj tobj(tid);
  
  rt.begin_epoch(0);
  
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], tid);
  }

  e.print();

  std::cout << "[===============================================================]" << std::endl;
  
  if (process) {
    e.process_global_buckets(tid);
    std::cout << "Num receives : " << rt.get_num_receives() << ", wis size : " << wis.size()
              << ", num buckets (syncs) : " << rt.get_num_syncs()
              << std::endl;

    assert(_numsyncs == rt.get_num_syncs());
    assert((rt.get_num_receives() == wis.size()) && "Number of receives is not equal to sends");
  }
}



template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm_multiple_threads(T1 _t1,
                                T2 _t2,
                                T3 _t3,
                                T4 _t4,
                                int _threads) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_configs<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(_threads, 1);
  
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);


  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], (i%_threads));
  }

  e.print();
}



template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm_multiple_numa_domains(T1 _t1,
                                     T2 _t2,
                                     T3 _t3,
                                     T4 _t4,
                                     int _threads,
                                     int _numa) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_configs<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(_threads, _numa);
  
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);

  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], (i%_threads));
  }

  e.print();
}




template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm_multiple_numa_domains_multiple_threads(T1 _t1,
                                                      T2 _t2,
                                                      T3 _t3,
                                                      T4 _t4,
                                                      int _threads,
                                                      int _numa) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_configs<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);
  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(_threads, _numa);
  
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);


  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], i%_threads);
  }

  //e.print();
}


template<typename buckets, typename runtime>
class ampp_thread_executor {

private:
  int nthreads;
  buckets& all_bkts;
  runtime& rt;
  bool process;  
  
public:
  ampp_thread_executor(int _n,
                       buckets& _all,
                       runtime& _rt,
                       bool _p) : nthreads(_n),
                                  all_bkts(_all),
                                  rt(_rt),
                                  process(_p){}
  
  void operator()(int tid,
                  std::vector<work_item_t>& wis) {
    // amplusplus::detail::push_thread_id_obj tobj(tid);
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
    
      for (typename std::vector<work_item_t>::size_type i = tid ;
           i < wis.size(); i+= nthreads) {
        work_item_t& wi = wis[i];
        //      debug("Pushing workitems in tid : ", tid);
        all_bkts.push(wi, tid);      
      }

      if (process) {
        //      debug("Processing in tid : ", tid);
        all_bkts.process_global_buckets(tid);
      }
    }
  }  
};


template<typename ordered_config_t>
void mt_general_test_eagm(ordered_config_t& config,
                          int _threads,
                          int _numa,
                          amplusplus::transport& trans,
                          int _numsyncs,
                          bool process) {



  CHAOTIC_ORDERING_T ch;
  debug("Setting threads to : ", _threads);
  trans.set_nthreads(_threads);
  typedef test_runtime<work_item_t> runtime_t;
  runtime_t rt(_threads, _numa, trans);

  debug("Runtime initialized ...");
  typedef boost::graph::agm::eagm_traits<typename ordered_config_t::global_ordering_t,
                                         typename ordered_config_t::node_ordering_t,
                                         typename ordered_config_t::numa_ordering_t,
                                         typename ordered_config_t::thread_ordering_t,
                                         ordered_config_t,
                                         runtime_t> eagm_traits_t;

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename eagm_traits_t::root_level_bucket_trait> all_buckets_t;
  all_buckets_t e(config, rt);

  std::vector<work_item_t> wis;
  init_data(wis);

  ampp_thread_executor<all_buckets_t, runtime_t> te(_threads, e, rt, process);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[_threads - 1]);

  for (int i = 0; i < _threads - 1; ++i) {
    boost::thread thr(boost::ref(te), i+1, wis);
    threads[i].swap(thr);
  }
	  
  te(0, wis);
    
  for (int i = 0; i < (_threads - 1); ++i)
    threads[i].join();


  trans.set_nthreads(1);
  e.print();

  if (process) {
    std::cout << "Num receives : " << rt.get_num_receives() << ", wis size : " << wis.size()
              << ", num buckets (syncs) : " << rt.get_num_syncs()
              << std::endl;

    assert(_numsyncs == rt.get_num_syncs());
    assert((rt.get_num_receives() == wis.size()) && "Number of receives is not equal to sends");
  }
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void mt_test_eagm(T1 _t1,
                  T2 _t2,
                  T3 _t3,
                  T4 _t4,
                  int _threads,
                  int _numa,
                  amplusplus::transport& trans,
                  int _numsyncs,
                  bool process = false) {


  std::string name = "[Global:";/* + std::string(_t1.name()) + 
    + ", Node:" + std::string(_t2.name()) + ", Numa:" + std::string(_t3.name())
    + ", Thread:" + std::string(_t4.name()) + "]";*/
  
  std::cout << "Running : "
            << name << ", Threads : " << _threads
            << ", Numa domains : " << _numa
            << ", Process : " << process 
            << std::endl;
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_configs<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);
  
  mt_general_test_eagm(config, _threads, _numa, trans, _numsyncs, process);
}


template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void mt_test_eagm_with_node_pq(T1 _t1,
                               T2 _t2,
                               T3 _t3,
                               T4 _t4,
                               int _threads,
                               int _numa,
                               amplusplus::transport& trans,
                               int _numsyncs,
                               bool process = false) {


  std::string name = "[Global:";/* + std::string(_t1.name()) + 
    + ", Node:" + std::string(_t2.name()) + ", Numa:" + std::string(_t3.name())
    + ", Thread:" + std::string(_t4.name()) + "]";*/
  
  std::cout << "Running With Node PQ : "
            << name << ", Threads : " << _threads
            << ", Numa domains : " << _numa
            << ", Process : " << process 
            << std::endl;
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4,
                                          boost::graph::agm::pq_container> ordered_config_t;

  ordered_config_t config = create_eagm_config_w_node_pq<T1,
                                                         T2,
                                                         T3,
                                                         T4>(_t1,
                                                             _t2,
                                                             _t3,
                                                             _t4);
  
  mt_general_test_eagm(config, _threads, _numa, trans, _numsyncs, process);
}


template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void mt_test_eagm_with_numa_pq(T1 _t1,
                               T2 _t2,
                               T3 _t3,
                               T4 _t4,
                               int _threads,
                               int _numa,
                               amplusplus::transport& trans,
                               int _numsyncs,
                               bool process = false) {


  std::string name = "[Global:";/* + std::string(_t1.name()) + 
    + ", Node:" + std::string(_t2.name()) + ", Numa:" + std::string(_t3.name())
    + ", Thread:" + std::string(_t4.name()) + "]";*/
  
  std::cout << "Running With Node PQ : "
            << name << ", Threads : " << _threads
            << ", Numa domains : " << _numa
            << ", Process : " << process 
            << std::endl;
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4,
                                          boost::graph::agm::buffer_container,
                                          boost::graph::agm::pq_container> ordered_config_t;

  ordered_config_t config = create_eagm_config_w_numa_pq<T1,
                                                         T2,
                                                         T3,
                                                         T4>(_t1,
                                                             _t2,
                                                             _t3,
                                                             _t4);
  
  mt_general_test_eagm(config, _threads, _numa, trans, _numsyncs, process);
}


int main(int argc, char* argv[]) {

  debug("Creating environment ...");
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  debug("Creating transport ...");
  amplusplus::transport trans = env.create_transport();
  debug("Transport created ...");
    
  CHAOTIC_ORDERING_T ch;
  DIJKSTRA_ORDERING_T dj;
  DIJKSTRA_ORDERING_STD_PQ_T dj_std_pq;
  DELTA_ORDERING_T delta(10);
  DELTA_ORDERING_T delta1(20);
  DELTA_ORDERING_T delta2(5); 

  std::cout << "Starting tests ..." << std::endl;
  
  test_eagm(ch, ch, ch, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(delta, ch, ch, ch, trans);  
  std::cout << "======================================================" << std::endl;
  test_eagm(delta, ch, ch, dj, trans);
  std::cout << "======================================================" << std::endl;  
  test_eagm(delta, delta2, ch, ch, trans);
  std::cout << "======================================================" << std::endl;    
  test_eagm(delta, delta2, ch, dj, trans);
  std::cout << "======================================================" << std::endl;    
  test_eagm(delta, ch, delta2, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(delta, delta2, dj, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(ch, delta, delta2, dj_std_pq, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(delta1, delta, delta2, dj_std_pq, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(ch, ch, ch, dj_std_pq, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(delta, ch, ch, dj_std_pq, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(ch, delta, ch, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(ch, dj, ch, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm_with_node_pq(ch, dj_std_pq, ch, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(ch, delta, ch, dj, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm(ch, ch, dj, ch, trans);
  std::cout << "======================================================" << std::endl;
  test_eagm_with_numa_pq(ch, ch, dj_std_pq, ch, trans);
  std::cout << "======================================================" << std::endl;
  

  int threads = 4;
  
  //test_eagm_multiple_threads(ch, ch, ch, ch, threads);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_threads(ch, ch, ch, dj_std_pq, threads);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_threads(delta, ch, ch, dj_std_pq, threads);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_threads(ch, delta, ch, ch, threads);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_threads(ch, delta, ch, dj, threads);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_threads(ch, ch, dj, ch, threads);
  //std::cout << "======================================================" << std::endl;

  int numa_domains = 2;

  //test_eagm_multiple_numa_domains(ch, ch, ch, ch, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;  
  //test_eagm_multiple_numa_domains(ch, ch, dj, ch, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;  
  //test_eagm_multiple_numa_domains(delta, ch, ch, dj_std_pq, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_numa_domains(ch, delta, ch, ch, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;
  //test_eagm_multiple_numa_domains(ch, delta, dj, ch, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;  
  //test_eagm_multiple_numa_domains(ch, ch, delta, dj, threads, numa_domains);

  numa_domains = 2;
  //std::cout << "======================================================" << std::endl;  
  //mt_test_eagm(ch, ch, delta, dj, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;  
  //mt_test_eagm(ch, ch, ch, ch, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;
  //mt_test_eagm(ch, ch, ch, dj, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;  
  //mt_test_eagm(delta, ch, ch, dj, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;
  //mt_test_eagm(ch, delta, ch, ch, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;
  //mt_test_eagm(ch, delta, ch, dj, threads, numa_domains);
  //std::cout << "======================================================" << std::endl;  
  //mt_test_eagm(ch, ch, dj, ch, threads, numa_domains);

  // With processing
  //test_eagm(ch, ch, ch, ch, trans, 1, true);
  //std::cout << "======================================================" << std::endl;
  //test_eagm(delta, ch, ch, ch, trans,5, true);  
  //std::cout << "======================================================" << std::endl;
  //test_eagm(delta, ch, ch, dj, trans, 5, true);
  //std::cout << "======================================================" << std::endl;  
  //test_eagm(delta, delta2, ch, ch, trans, 5, true);
  //std::cout << "======================================================" << std::endl;    
  //test_eagm(delta, delta2, ch, dj, trans, 5, true);
  // std::cout << "======================================================" << std::endl;    
  // test_eagm(delta, ch, delta2, ch, trans, 5, true);
  //    std::cout << "======================================================" << std::endl;
  //test_eagm(delta, delta2, dj, ch, trans, 5, true);
  //std::cout << "======================================================" << std::endl;
  //test_eagm(ch, delta, delta2, dj_std_pq, trans, 1, true);
  //std::cout << "======================================================" << std::endl;
  //test_eagm(delta1, delta, delta2, dj_std_pq, trans, 3, true);
  // std::cout << "======================================================" << std::endl;
  //test_eagm(ch, ch, ch, dj_std_pq, trans, 1, true);
  //    std::cout << "======================================================" << std::endl;
  //test_eagm(delta, ch, ch, dj_std_pq, trans, 5, true);
  //    std::cout << "======================================================" << std::endl;
  //test_eagm(ch, delta, ch, ch, trans, 1, true);
  //    std::cout << "======================================================" << std::endl;
  //  test_eagm(ch, dj, ch, ch, trans, 1, true);
  //    std::cout << "======================================================" << std::endl;
  //  test_eagm_with_node_pq(ch, dj_std_pq, ch, ch, trans, 1, true);
  //    std::cout << "======================================================" << std::endl;
  //  test_eagm(ch, delta, ch, dj, trans, 1, true);
  //     std::cout << "======================================================" << std::endl;
  //   test_eagm(ch, ch, dj, ch, trans, 1, true);
  //    std::cout << "======================================================" << std::endl;
  //  test_eagm_with_numa_pq(ch, ch, dj_std_pq, ch, trans, 1, true);
  //  std::cout << "======================================================" << std::endl;

  threads = 8;
  numa_domains = 1;

  int iterations = 3;
  
  for (int i=0; i < iterations; ++i) {
    for (int i=0; i < 4; ++i) {

      threads = std::pow(2, i);
      numa_domains = 1;
    
      mt_test_eagm(ch, ch, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, delta, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, ch, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta, ch, ch, dj, threads, numa_domains, trans, 5, true);
      mt_test_eagm(ch, delta, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, delta, ch, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, dj, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta, ch, ch, ch, threads, numa_domains, trans,5, true);
      mt_test_eagm(delta, delta2, ch, ch, threads, numa_domains, trans, 5, true);
      mt_test_eagm(delta, delta2, ch, dj, threads, numa_domains, trans, 5, true);
      mt_test_eagm(delta, ch, delta2, ch, threads, numa_domains, trans, 5, true);
      mt_test_eagm(delta, delta2, dj, ch, threads, numa_domains, trans, 5, true);
      mt_test_eagm(ch, delta, delta2, dj_std_pq, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta1, delta, delta2, dj_std_pq, threads, numa_domains, trans, 3, true);
      mt_test_eagm(ch, ch, ch, dj_std_pq, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta, ch, ch, dj_std_pq, threads, numa_domains, trans, 5, true);
      mt_test_eagm(ch, delta, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, dj, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, delta, ch, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, dj, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm_with_node_pq(ch, dj_std_pq, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm_with_numa_pq(ch, ch, dj_std_pq, ch, threads, numa_domains, trans, 1, true);
    }

    for (int i=1; i < 4; ++i) {

      threads = std::pow(2, i);
      numa_domains = 2;
    
      mt_test_eagm(ch, ch, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, delta, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, ch, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta, ch, ch, dj, threads, numa_domains, trans, 5, true);
      mt_test_eagm(ch, delta, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, delta, ch, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, dj, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta, ch, ch, ch, threads, numa_domains, trans,5, true);
      mt_test_eagm(delta, delta2, ch, ch, threads, numa_domains, trans, 5, true);
      mt_test_eagm(delta, delta2, ch, dj, threads, numa_domains, trans, 5, true);
      mt_test_eagm(delta, ch, delta2, ch, threads, numa_domains, trans, 5, true);
      mt_test_eagm(delta, delta2, dj, ch, threads, numa_domains, trans, 5, true);
      mt_test_eagm(ch, delta, delta2, dj_std_pq, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta1, delta, delta2, dj_std_pq, threads, numa_domains, trans, 3, true);
      mt_test_eagm(ch, ch, ch, dj_std_pq, threads, numa_domains, trans, 1, true);
      mt_test_eagm(delta, ch, ch, dj_std_pq, threads, numa_domains, trans, 5, true);
      mt_test_eagm(ch, delta, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, dj, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, delta, ch, dj, threads, numa_domains, trans, 1, true);
      mt_test_eagm(ch, ch, dj, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm_with_node_pq(ch, dj_std_pq, ch, ch, threads, numa_domains, trans, 1, true);
      mt_test_eagm_with_numa_pq(ch, ch, dj_std_pq, ch, threads, numa_domains, trans, 1, true);
    }
  }
    
}
