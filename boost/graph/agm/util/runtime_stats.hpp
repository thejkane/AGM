#ifndef __AGM_RUNTIME_STATS__
#define __AGM_RUNTIME_STATS__
class runtime_stats {
private:
  int n_all_to_alls;
  uint64_t n_epochs;  
  std::vector<uint64_t> n_sends;
  std::vector<uint64_t> n_receives;
  int threads;

public:
  runtime_stats(int _threads): threads(_threads),
                               n_all_to_alls(0),
                               n_epochs(0){
    //    std::cout << "threads in constructor : " << threads << std::endl;
    n_sends.resize(threads, 0);
    n_receives.resize(threads, 0);
  }

  void increment_all_to_alls() {
#ifdef PBGL2_WORK_STATS
    n_all_to_alls++;
#endif
  }

  void increment_epochs(int tid) {
#ifdef PBGL2_WORK_STATS
    if (tid == 0)
      n_epochs++;
#endif
  }

  void increment_sends(int tid) {
    //std::cout << "tid in increment send : " << tid << std::endl;    
#ifdef PBGL2_WORK_STATS
    n_sends[tid]++;
#endif
  }
  
  void increment_receives(int tid) {
#ifdef PBGL2_WORK_STATS
    n_receives[tid]++;
#endif
  }  
};
#endif
