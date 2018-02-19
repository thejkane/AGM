#include "../common/utils.hpp"
#include <boost/graph/agm/util/eagm_buckets.hpp>
#include <boost/tuple/tuple.hpp>

#include <limits.h>
// These tests are for dynamic configuration

struct chaotic {
public:
  template <typename T>
  bool operator()(T i, T j) {
    return false;
  }
};

template<int index>
struct dijkstra {
public:
  template <typename T>
  bool operator()(T i, T j) {
    return (std::get<index>(i) < std::get<index>(j));
  }
};

template<int index, typename delta_t>
struct delta_ord {

private:
  delta_t delta;
  
public:
  delta_ord(delta_t _d) : delta(_d) {}

  delta_ord(const delta_ord& _do) : delta(_do.delta) {}
  
  template <typename T>
  bool operator()(T i, T j) {
    return ((std::get<index>(i)/delta) < (std::get<index>(j)/delta));
  }
};


struct chaotic_ordering_gen {
public:
  typedef chaotic strict_ordering;
};

//typedef boost::tuple<uint64_t, uint64_t> pwork_item_t;
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

template<typename ordering>
void execute(int argc, char* argv[],
             ordering _ord,
             uint64_t expected) {

  std::vector<work_item_t> wis;
  init_data(wis);
    
  /*  typedef boost::graph::agm::buckets<work_item_t, ordering> buckets_t;
  buckets_t all_buckets(_ord);

  for (int i=0; i < wis.size(); i++) {
    all_buckets.push(wis[i]);
  }

  assert(all_buckets.get_buckets_created() == expected);
  all_buckets.print_buckets();*/
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
      all_bkts.push(wi, numa, tid);      
    }  
  }  
};

template<typename ordering>
void execute_mt(int argc, char* argv[],
                ordering _ord,
                int nthreads,
                uint64_t expected) {

  std::vector<work_item_t> wis;
  init_data(wis);

  
  /*typedef boost::graph::agm::buckets<work_item_t, ordering> buckets_t;
  buckets_t all_buckets(_ord);

  thread_executor<buckets_t> te(nthreads, all_buckets);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
  for (int i = 0; i < nthreads - 1; ++i) {
    boost::thread thr(boost::ref(te), i + 1, wis);
    threads[i].swap(thr);
  }
	  
  te(0, wis);
    
  for (int i = 0; i < (nthreads - 1); ++i)
    threads[i].join();

  std::cout << "All buckets created = " << all_buckets.get_buckets_created()
            << std::endl;
  
  assert(all_buckets.get_buckets_created() == expected);
  all_buckets.print_buckets();*/
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
boost::graph::agm::eagm_configs<T1, T2, T3, T4> create_eagm_config(T1 _t1,
                                                                   T2 _t2,
                                                                   T3 _t3,
                                                                   T4 _t4) {
  typedef boost::graph::agm::eagm_configs<T1, T2, T3, T4> eagm_config_t;
  eagm_config_t config(_t1,
                       _t2,
                       _t3,
                       _t4);

  config.ranks = 1;
  config.threads = 1;
  config.numa_domains = 1;

  return config;
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_optimizations(T1 _t1,
                        T2 _t2,
                        T3 _t3,
                        T4 _t4) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_config<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  config.optimize();
  config.print();                             
}


template<typename T1,
         typename T2,
         typename T3,
         typename T4>
void test_eagm(T1 _t1,
               T2 _t2,
               T3 _t3,
               T4 _t4) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_config<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  config.optimize();
  config.print();

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename ordered_config_t::global_ordering_t,
                                          ordered_config_t> eagms;

  eagms e(config.global_ord,
          config,
          boost::graph::agm::spatial_level::global);

  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], 0, 0);
  }

  e.print();
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
  ordered_config_t config = create_eagm_config<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  config.threads = _threads;
  config.optimize();
  config.print();

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename ordered_config_t::global_ordering_t,
                                          ordered_config_t> eagms;

  eagms e(config.global_ord,
          config,
          boost::graph::agm::spatial_level::global);

  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], 0, (i%config.threads));
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
                                int _numa) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_config<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  config.numa_domains = _numa;
  config.optimize();
  config.print();

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename ordered_config_t::global_ordering_t,
                                          ordered_config_t> eagms;

  eagms e(config.global_ord,
          config,
          boost::graph::agm::spatial_level::global);

  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], (i%config.numa_domains), 0);
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
  ordered_config_t config = create_eagm_config<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);
  //configs
  config.numa_domains = _numa;
  config.threads = _threads;
  config.optimize();
  config.print();

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename ordered_config_t::global_ordering_t,
                                          ordered_config_t> eagms;

  eagms e(config.global_ord,
          config,
          boost::graph::agm::spatial_level::global);

  std::vector<work_item_t> wis;
  init_data(wis);
    
  for (int i=0; i < wis.size(); i++) {
    e.push(wis[i], (i%config.numa_domains), i%config.threads);
  }

  e.print();
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
                  int _numa) {
  
  typedef boost::graph::agm::eagm_configs<T1,
                                          T2,
                                          T3,
                                          T4> ordered_config_t;

  CHAOTIC_ORDERING_T ch;
  ordered_config_t config = create_eagm_config<T1,
                                               T2,
                                               T3,
                                               T4>(_t1,
                                                   _t2,
                                                   _t3,
                                                   _t4);

  config.numa_domains = _numa;
  config.threads = _threads;
  config.optimize();
  config.print();

  typedef boost::graph::agm::eagm_buckets<work_item_t,
                                          typename ordered_config_t::global_ordering_t,
                                          ordered_config_t> eagms;

  eagms e(config.global_ord,
          config,
          boost::graph::agm::spatial_level::global);

  std::vector<work_item_t> wis;
  init_data(wis);

  thread_executor<eagms> te(config.threads, e);
  
  boost::scoped_array<boost::thread> threads(new boost::thread[config.threads - 1]);
  for (int i = 0; i < config.threads - 1; ++i) {
    boost::thread thr(boost::ref(te), (i%config.numa_domains), i+1, wis);
    threads[i].swap(thr);
  }
	  
  te(0, 0, wis);
    
  for (int i = 0; i < (config.threads - 1); ++i)
    threads[i].join();

  e.print();
}


int main(int argc, char* argv[]) {
  
  CHAOTIC_ORDERING_T ch;
  DIJKSTRA_ORDERING_T dj;
  DELTA_ORDERING_T delta(10); 
  /*  
      test_optimizations(ch, ch, ch, ch);

      std::cout << "======================================================" << std::endl;
  
      test_optimizations(ch, ch, ch, dj);

      std::cout << "======================================================" << std::endl;
  
      test_optimizations(delta, ch, ch, dj);

      std::cout << "======================================================" << std::endl;
      test_optimizations(ch, delta, ch, dj);

      std::cout << "======================================================" << std::endl;
      test_optimizations(ch, delta, dj, dj);

      std::cout << "======================================================" << std::endl;
      test_optimizations(ch, ch, dj, ch);*/


  /*test_eagm(ch, ch, ch, ch);

  std::cout << "======================================================" << std::endl;
  
  test_eagm(ch, ch, ch, dj);

  std::cout << "======================================================" << std::endl;
  
  test_eagm(delta, ch, ch, dj);

  std::cout << "======================================================" << std::endl;

  test_eagm(ch, delta, ch, ch);

  std::cout << "======================================================" << std::endl;

  test_eagm(ch, delta, ch, dj);

  std::cout << "======================================================" << std::endl;
  
  test_eagm(ch, ch, dj, ch);*/


  int threads = 4;
  int numa_domains = 2;
  
  //test_eagm_multiple_threads(ch, ch, ch, ch, threads);

  //std::cout << "======================================================" << std::endl;
  
  //test_eagm_multiple_threads(ch, ch, ch, dj, threads);

  //std::cout << "======================================================" << std::endl;
  
  //test_eagm_multiple_threads(delta, ch, ch, dj, threads);

  //std::cout << "======================================================" << std::endl;

  //test_eagm_multiple_threads(ch, delta, ch, ch, threads);

  //std::cout << "======================================================" << std::endl;

  //test_eagm_multiple_threads(ch, delta, ch, dj, threads);

  //std::cout << "======================================================" << std::endl;
  
  //test_eagm_multiple_threads(ch, ch, dj, ch, threads);





  //  test_eagm_multiple_numa_domains(ch, ch, ch, ch, numa_domains);

  //std::cout << "======================================================" << std::endl;
  
  //test_eagm_multiple_numa_domains(ch, ch, dj, ch, numa_domains);

  //std::cout << "======================================================" << std::endl;
  
  //test_eagm_multiple_numa_domains(delta, ch, ch, dj, numa_domains);

  //std::cout << "======================================================" << std::endl;

  //test_eagm_multiple_numa_domains(ch, delta, ch, ch, numa_domains);

  //std::cout << "======================================================" << std::endl;

  //test_eagm_multiple_numa_domains(ch, delta, dj, ch, numa_domains);

  //std::cout << "======================================================" << std::endl;
  
  //test_eagm_multiple_numa_domains(ch, ch, dj, ch, numa_domains);

  //test_eagm_multiple_numa_domains_multiple_threads(ch, ch, delta, dj, threads, numa_domains);

  mt_test_eagm(ch, ch, delta, dj, threads, numa_domains);

  mt_test_eagm(ch, ch, ch, ch, threads, numa_domains);

  std::cout << "======================================================" << std::endl;
  
  mt_test_eagm(ch, ch, ch, dj, threads, numa_domains);

  std::cout << "======================================================" << std::endl;
  
  mt_test_eagm(delta, ch, ch, dj, threads, numa_domains);

  std::cout << "======================================================" << std::endl;

  mt_test_eagm(ch, delta, ch, ch, threads, numa_domains);

  std::cout << "======================================================" << std::endl;

  mt_test_eagm(ch, delta, ch, dj, threads, numa_domains);

  std::cout << "======================================================" << std::endl;
  
  mt_test_eagm(ch, ch, dj, ch, threads, numa_domains);

}
