#ifndef __AGM_BUCKET_STATS__
#define __AGM_BUCKET_STATS__
#include <atomic>

class bucket_stats {
private:
  std::atomic<std::uint64_t> buckets_created;
  std::vector<std::uint64_t> total_pushes;

public:
  bucket_stats(int nthreads):buckets_created(0),
                             total_pushes(nthreads, 0){}
  void increment_buckets() {
    buckets_created++;
  }

  void increment_pushes(int tid) {
    ++total_pushes[tid];
  }

  void print() {
    std::uint64_t sum = std::accumulate(total_pushes.begin(), total_pushes.end(), 0);

    std::uint64_t global_sum = 0;
    MPI_Allreduce(&sum, &global_sum, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    info("Buckets created : ", buckets_created.load());
    info("Total pushes to the data-structure : ", global_sum);    
  }

  std::uint64_t get_number_of_bkst_created() {
    return (buckets_created.load());
  }
};

#endif
