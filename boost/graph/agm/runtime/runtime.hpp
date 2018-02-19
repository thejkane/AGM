#ifndef __EAGM_RUNTIME_INTERFACE__
#define __EAGM_RUNTIME_INTERFACE__

#include <functional>
#include <atomic>

namespace boost { namespace graph { namespace agm {

template<typename work_item>
class runtime {

public:

  void send(const work_item& wi, int tid) = 0;

  void initialize_per_thread(int tid) = 0;
  
  void register_receiver(std::function<void(const work_item&, int)> funcrev) = 0;

  int find_numa_node(int tid) = 0;

  int get_nnuma_nodes() = 0;
  
  int get_nthreads() = 0;

  int get_nranks() = 0;

  void do_all_gather(void* _ptosend,
                     int sendcnt,
                     void* _precv,
                     int recvcnt) = 0;

  void increase_activity_count(int tid, uint64_t v) = 0;
  
  void decrease_activity_count(int tid, uint64_t v) = 0;

  uint64_t get_activity_count() = 0;
  
  int get_nthreads_for_numa_domain(int tid) = 0;

  int get_thread_index_in_numa_domain(int tid) = 0;

  void wait_for_numa_domain_threads_to_reach_here(int tid) = 0;

  bool is_main_thread_in_numa_domain(int tid) = 0;
  
  bool begin_epoch(int tid) = 0;
  
  bool is_epoch_completed(int tid) = 0;

  void end_epoch() = 0;

  unsigned long end_epoch_with_value(const unsigned long& _read_value) = 0;

  void pull_work(int tid) = 0;
  
  void synchronize() = 0;

  void wait_for_threads_to_reach_here(int tid) = 0;

  bool is_main_thread(int tid) = 0;
};


class runtime_base {
protected:
  std::atomic<uint64_t> active_count;
  int threads;
  int numa_domains;
  int ranks;  

public:
  runtime_base(int _threads, int _numa, int _ranks):active_count(0),
                                                    threads(_threads),
                                                    numa_domains(_numa),
                                                    ranks(_ranks){}

  inline void increase_activity_count(int tid, uint64_t v) {
    active_count.fetch_add(v);
  }
  
  inline void decrease_activity_count(int tid, uint64_t v) {
    active_count.fetch_sub(v);
  }

  inline uint64_t get_activity_count() {
    return active_count.load();
  }

  int get_nnuma_nodes() {
    return numa_domains;
  }
  
  inline int get_nthreads() {
    return threads;
  }  

  inline int get_nranks() {
    return ranks;
  }
};

}}}
#endif
