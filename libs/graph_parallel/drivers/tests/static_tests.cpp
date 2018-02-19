#include <iostream>
#include <tuple>
#include <boost/tuple/tuple.hpp>
#include <limits.h>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/util/eagm_buckets.hpp>
#include <boost/graph/agm/model/bucket_handler.hpp>
#include <boost/graph/agm/util/eagm_config.hpp>
#include <boost/thread/locks.hpp>

template<typename work_item>
class test_runtime {
private:
  int threads;  
  int numa_domains;
  
public:

  typedef work_item work_item_t;
  
  test_runtime(int _nthreads,
               int _nnuma): threads(_nthreads),
                            numa_domains(_nnuma){}
  
  void send(work_item& wi) {}

  void initialize_per_thread(int tid){}
  
  void register_receiver(std::function<void(const work_item&, int)> funcrev) {
    // set message handler
  }

  int find_numa_node(int tid) {
    return (tid % numa_domains);
  }

  inline int get_nthreads() {
    return threads;
  }

  int get_nnuma_nodes() {
    return 1;
  }
  
  void increase_activity_count(int) {}
  void decrease_activity_count(int) {}
  
  void increase_mt_activity_count(int){}
  void decrease_mt_activity_count(int){}

  bool begin_epoch(){ return true; }
  bool is_epoch_completed() { return true; }

  void pull_work() {}
  
  void synchronize(){}

  void synchronize_threads(){}

  void synchronize_numa_domain_threads(int tid){}    
};


class processing_fn {
  template<typename work_item,
           typename buckets>
  void operator()(const work_item& wi,
                  int tid,
                  buckets _bkts) {
  }
};

int main() {

  using namespace boost::graph::agm;
  
  std::cout << "Starting ..." << std::endl;
  typedef std::tuple<int, int> work_item;
  work_item w(4,5);

  typedef test_runtime<work_item> runtime_t;
  runtime_t rt(1, 1);
  bucket_handler<CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 runtime_t,
                 processing_fn> chaotic;
  chaotic.print();

  DELTA_ORDERING_T delta(5);
  CHAOTIC_ORDERING_T ch;
  DIJKSTRA_ORDERING_T dj;

  processing_fn pf;

  typedef eagm_configs<DELTA_ORDERING_T,
               CHAOTIC_ORDERING_T,
               DELTA_ORDERING_T,
               DIJKSTRA_ORDERING_T> eagm_config_t;
    
  eagm_config_t eagmconfig(delta,
                           ch,
                           delta,
                           dj);
  
  bucket_handler<eagm_config_t::global_ordering_t,
                 eagm_config_t::node_ordering_t,
                 eagm_config_t::numa_ordering_t,
                 eagm_config_t::thread_ordering_t,
                 runtime_t,
                 processing_fn> general(pf, eagmconfig,
                                        rt);
  general.print();
  general.initialize();
  general.run(0);

  
  bucket_handler<DIJKSTRA_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 runtime_t,
                 processing_fn> djagm;

  djagm.print();

  bucket_handler<CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 CHAOTIC_ORDERING_T,
                 DIJKSTRA_ORDERING_T,
                 runtime_t,
                 processing_fn> dcagm;

  dcagm.print();
      

  {
    typedef eagm_configs<CHAOTIC_ORDERING_T,
                         CHAOTIC_ORDERING_T,
                         CHAOTIC_ORDERING_T,
                         DIJKSTRA_ORDERING_T> a_eagm_config_t;
    
    a_eagm_config_t conf(ch,
                         ch,
                         ch,
                         dj);

    // Global-->Thread-->Nil
    typedef boost::graph::agm::eagm_traits<CHAOTIC_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           DIJKSTRA_ORDERING_T,
                                           a_eagm_config_t,
                                           runtime_t> eagm_traits_t;

    typedef boost::graph::agm::eagm_buckets<work_item,
                                            eagm_traits_t::root_level_bucket_trait> all_buckets_t;
    all_buckets_t buckets(conf, rt);
    buckets.push(w, 0);
  }

  {

    typedef eagm_configs<CHAOTIC_ORDERING_T,
                         DIJKSTRA_ORDERING_T,
                         CHAOTIC_ORDERING_T,
                         CHAOTIC_ORDERING_T> a_eagm_config_t;
    
    a_eagm_config_t conf(ch,
                         dj,
                         ch,
                         ch);

    // Global-->Node-->Nil
    typedef boost::graph::agm::eagm_traits<CHAOTIC_ORDERING_T,
                                           DIJKSTRA_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           a_eagm_config_t,
                                           runtime_t> eagm_traits_t;

    typedef boost::graph::agm::eagm_buckets<work_item,
                                            eagm_traits_t::root_level_bucket_trait> all_buckets_t;
    all_buckets_t::spatial_t spatiallevel;
    all_buckets_t buckets(conf, rt);
    buckets.push(w, 0);
  }
  
  {
    typedef eagm_configs<CHAOTIC_ORDERING_T,
                         CHAOTIC_ORDERING_T,
                         DIJKSTRA_ORDERING_T,
                         CHAOTIC_ORDERING_T> a_eagm_config_t;
    
    a_eagm_config_t conf(ch,
                         ch,
                         dj,
                         ch);

    // Global-->Numa-->Nil
    typedef boost::graph::agm::eagm_traits<CHAOTIC_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           DIJKSTRA_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           a_eagm_config_t,
                                           runtime_t> eagm_traits_t;

    typedef boost::graph::agm::eagm_buckets<work_item,
                                            eagm_traits_t::root_level_bucket_trait> all_buckets_t;

    all_buckets_t buckets(conf, rt);
    buckets.push(w, 0);
  }

  {
    typedef eagm_configs<CHAOTIC_ORDERING_T,
                         DELTA_ORDERING_T,
                         CHAOTIC_ORDERING_T,
                         DIJKSTRA_ORDERING_T> a_eagm_config_t;
    
    a_eagm_config_t conf(ch,
                         delta,
                         ch,
                         dj);

    typedef boost::graph::agm::eagm_traits<CHAOTIC_ORDERING_T,
                                           DELTA_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           DIJKSTRA_ORDERING_T,
                                           a_eagm_config_t,
                                           runtime_t> eagm_traits_t;

    typedef boost::graph::agm::eagm_buckets<work_item,
                                            eagm_traits_t::root_level_bucket_trait> all_buckets_t;
    all_buckets_t buckets(conf, rt);
    buckets.push(w, 0);
  }

  {
    typedef eagm_configs<CHAOTIC_ORDERING_T,
                         DELTA_ORDERING_T,
                         DELTA_ORDERING_T,
                         DIJKSTRA_ORDERING_T> a_eagm_config_t;
    
    a_eagm_config_t conf(ch,
                         delta,
                         delta,
                         dj);

    // 
    typedef boost::graph::agm::eagm_traits<CHAOTIC_ORDERING_T,
                                           DELTA_ORDERING_T,
                                           DELTA_ORDERING_T,
                                           DIJKSTRA_ORDERING_T,
                                           a_eagm_config_t,
                                           runtime_t> eagm_traits_t;

    typedef boost::graph::agm::eagm_buckets<work_item,
                                            eagm_traits_t::root_level_bucket_trait> all_buckets_t;
    all_buckets_t buckets(conf, rt);
    buckets.push(w, 0);

  }

  {
    typedef eagm_configs<CHAOTIC_ORDERING_T,
                         CHAOTIC_ORDERING_T,
                         DELTA_ORDERING_T,
                         DIJKSTRA_ORDERING_T> a_eagm_config_t;
    
    a_eagm_config_t conf(ch,
                         ch,
                         delta,
                         dj);

    
    // 
    typedef boost::graph::agm::eagm_traits<CHAOTIC_ORDERING_T,
                                           CHAOTIC_ORDERING_T,
                                           DELTA_ORDERING_T,
                                           DIJKSTRA_ORDERING_T,
                                           a_eagm_config_t,
                                           runtime_t> eagm_traits_t;

    typedef boost::graph::agm::eagm_buckets<work_item,
                                            eagm_traits_t::root_level_bucket_trait> all_buckets_t;
    all_buckets_t buckets(conf, rt);
    buckets.push(w, 0);

  }
}
