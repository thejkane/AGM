#ifndef PBGL2_WORK_STATS
#define PBGL2_WORK_STATS

class iteration_stats {
public:
  uint64_t useful;
  uint64_t invalidated;
  uint64_t invalidated_cancels;  
  uint64_t rejected;
  uint64_t edges;
  double teps;

  iteration_stats():useful(0),
		    invalidated(0),
                    invalidated_cancels(0),
		    rejected(0),
		    edges(0),
		    teps((double)0){}

  iteration_stats(const iteration_stats &is) : useful(is.useful),
					       invalidated(is.invalidated),
                                               invalidated_cancels(is.invalidated_cancels),
					       rejected(is.rejected),
					       edges(is.edges),
					       teps(is.teps) {}
};

class agm_work_stats {
  typedef std::tuple<uint64_t, uint64_t, uint64_t , uint64_t, uint64_t> work_stats_t;


private:
  int threads;
  std::vector<uint64_t> useful;
  std::vector<uint64_t> invalidated;
  std::vector<uint64_t> invalidated_cancels;  
  std::vector<uint64_t> rejected;
  std::vector<uint64_t> edges_traversed;

  std::vector<iteration_stats> stats_per_iteration;
public:
  agm_work_stats(int _t) : threads(_t) {
    useful.resize(threads, 0);
    invalidated.resize(threads, 0);
    invalidated_cancels.resize(threads, 0);
    rejected.resize(threads, 0);
    edges_traversed.resize(threads, 0);
  }

  void increment_useful(int tid) {
#ifdef PBGL2_PRINT_WORK_STATS
    ++useful[tid];
#endif
  }

  void increment_invalidated(int tid) {
#ifdef PBGL2_PRINT_WORK_STATS
    ++invalidated[tid];
#endif
  }

  void increment_invalidated_cancels(int tid) {
#ifdef PBGL2_PRINT_WORK_STATS
    ++invalidated_cancels[tid];
#endif    
  }

  void increment_rejected(int tid) {
#ifdef PBGL2_PRINT_WORK_STATS
    ++rejected[tid];
#endif
  }

  void increment_edges(int tid) {
#ifdef PBGL2_PRINT_WORK_STATS
    ++edges_traversed[tid];
#endif
  }

  void reset() {
    for(unsigned int i = 0; i < useful.size(); ++i) {
      useful[i] = 0;
      invalidated[i] = 0;
      invalidated_cancels[i] = 0;
      rejected[i] = 0;
      edges_traversed[i] = 0;
    }   
  }

  void reduce_stats(time_type _executiont) {
#ifdef PBGL2_PRINT_WORK_STATS
    // first reduce locally
    work_stats_t stats{ (uint64_t)0, (uint64_t)0, (uint64_t)0, (uint64_t)0, (uint64_t)0 };
    for(unsigned int i = 0; i < useful.size(); ++i) {
      std::get<0>(stats) += useful[i];
      std::get<1>(stats) += invalidated[i];
      std::get<2>(stats) += invalidated_cancels[i];      
      std::get<3>(stats) += rejected[i];
      std::get<4>(stats) += edges_traversed[i];
    }
    
    // reduce globally
    work_stats_t temp_stats;
    MPI_Allreduce(&std::get<0>(stats), &std::get<0>(temp_stats), 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&std::get<1>(stats), &std::get<1>(temp_stats), 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&std::get<2>(stats), &std::get<2>(temp_stats), 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);    
    MPI_Allreduce(&std::get<3>(stats), &std::get<3>(temp_stats), 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&std::get<4>(stats), &std::get<4>(temp_stats), 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    iteration_stats istats;
    istats.useful = std::get<0>(temp_stats);
    istats.invalidated = std::get<1>(temp_stats);
    istats.invalidated_cancels = std::get<2>(temp_stats);
    istats.rejected = std::get<3>(temp_stats);
    istats.edges = std::get<4>(temp_stats);
    istats.teps = (double)(std::get<4>(temp_stats) / _executiont);

    // accumulate into an iteration
    stats_per_iteration.push_back(istats);

    // reset counters
    reset();
#endif
  }

  void print_summary() {
#ifdef PBGL2_PRINT_WORK_STATS
    std::vector<double> allteps;
    //stats
    std::vector<uint64_t> alluseful;
    std::vector<uint64_t> allinvalidated;
    std::vector<uint64_t> allinvalidatedcancels;    
    std::vector<uint64_t> allrejected;
    std::vector<uint64_t> alledges;

    for(int i=0; i < stats_per_iteration.size(); ++i) {
      alluseful.push_back(stats_per_iteration[i].useful);
      allinvalidated.push_back(stats_per_iteration[i].invalidated);
      allinvalidatedcancels.push_back(stats_per_iteration[i].invalidated_cancels);      
      allrejected.push_back(stats_per_iteration[i].rejected);
      alledges.push_back(stats_per_iteration[i].edges);
      allteps.push_back(stats_per_iteration[i].teps);
    }

    // TEPS
    double mean_teps = 0;
    double min_teps = 0;
    double q1_teps = 0;
    double median_teps = 0;
    double q3_teps = 0;
    double max_teps = 0;
    double stddev_teps = 0;

    time_statistics(allteps,
		    mean_teps,
		    min_teps,
		    q1_teps,
		    median_teps,
		    q3_teps,
		    max_teps,
		    stddev_teps);

    std::cout << "[TEPS SUMMARY] MEAN : " << mean_teps
	      << ", MIN : " << min_teps
	      << ", Q1 : " << q1_teps
	      << ", MEDIAN : " << median_teps
	      << ", Q3 : " << q3_teps
	      << ", MAX : " << max_teps
	      << ", STDDEV : " << stddev_teps
	      << std::endl;

    // USEFUL
    uint64_t mean_useful = 0;
    uint64_t min_useful = 0;
    uint64_t q1_useful = 0;
    uint64_t median_useful = 0;
    uint64_t q3_useful = 0;
    uint64_t max_useful = 0;
    uint64_t stddev_useful = 0;

    time_statistics(alluseful,
		    mean_useful,
		    min_useful,
		    q1_useful,
		    median_useful,
		    q3_useful,
		    max_useful,
		    stddev_useful);

    std::cout << "[USEFUL SUMMARY] MEAN : " << mean_useful
	      << ", MIN : " << min_useful
	      << ", Q1 : " << q1_useful
	      << ", MEDIAN : " << median_useful
	      << ", Q3 : " << q3_useful
	      << ", MAX : " << max_useful
	      << ", STDDEV : " << stddev_useful
	      << std::endl;

    // INVALIDATED
    uint64_t mean_invalidated = 0;
    uint64_t min_invalidated = 0;
    uint64_t q1_invalidated = 0;
    uint64_t median_invalidated = 0;
    uint64_t q3_invalidated = 0;
    uint64_t max_invalidated = 0;
    uint64_t stddev_invalidated = 0;

    time_statistics(allinvalidated,
		    mean_invalidated,
		    min_invalidated,
		    q1_invalidated,
		    median_invalidated,
		    q3_invalidated,
		    max_invalidated,
		    stddev_invalidated);

    std::cout << "[INVALIDATED SUMMARY] MEAN : " << mean_invalidated
	      << ", MIN : " << min_invalidated
	      << ", Q1 : " << q1_invalidated
	      << ", MEDIAN : " << median_invalidated
	      << ", Q3 : " << q3_invalidated
	      << ", MAX : " << max_invalidated
	      << ", STDDEV : " << stddev_invalidated
	      << std::endl;

    // INVALIDATED CANCELS
    uint64_t mean_invalidatedcancels = 0;
    uint64_t min_invalidatedcancels = 0;
    uint64_t q1_invalidatedcancels = 0;
    uint64_t median_invalidatedcancels = 0;
    uint64_t q3_invalidatedcancels = 0;
    uint64_t max_invalidatedcancels = 0;
    uint64_t stddev_invalidatedcancels = 0;

    time_statistics(allinvalidatedcancels,
		    mean_invalidatedcancels,
		    min_invalidatedcancels,
		    q1_invalidatedcancels,
		    median_invalidatedcancels,
		    q3_invalidatedcancels,
		    max_invalidatedcancels,
		    stddev_invalidatedcancels);

    std::cout << "[INVALIDATED CANCELS SUMMARY] MEAN : " << mean_invalidatedcancels
	      << ", MIN : " << min_invalidatedcancels
	      << ", Q1 : " << q1_invalidatedcancels
	      << ", MEDIAN : " << median_invalidatedcancels
	      << ", Q3 : " << q3_invalidatedcancels
	      << ", MAX : " << max_invalidatedcancels
	      << ", STDDEV : " << stddev_invalidatedcancels
	      << std::endl;
    
    
    // REJECTED
    uint64_t mean_rejected = 0;
    uint64_t min_rejected = 0;
    uint64_t q1_rejected = 0;
    uint64_t median_rejected = 0;
    uint64_t q3_rejected = 0;
    uint64_t max_rejected = 0;
    uint64_t stddev_rejected = 0;

    time_statistics(allrejected,
		    mean_rejected,
		    min_rejected,
		    q1_rejected,
		    median_rejected,
		    q3_rejected,
		    max_rejected,
		    stddev_rejected);

    std::cout << "[REJECTED SUMMARY] MEAN : " << mean_rejected
	      << ", MIN : " << min_rejected
	      << ", Q1 : " << q1_rejected
	      << ", MEDIAN : " << median_rejected
	      << ", Q3 : " << q3_rejected
	      << ", MAX : " << max_rejected
	      << ", STDDEV : " << stddev_rejected
	      << std::endl;

    // EDGES
    uint64_t mean_edges = 0;
    uint64_t min_edges = 0;
    uint64_t q1_edges = 0;
    uint64_t median_edges = 0;
    uint64_t q3_edges = 0;
    uint64_t max_edges = 0;
    uint64_t stddev_edges = 0;

    time_statistics(alledges,
		    mean_edges,
		    min_edges,
		    q1_edges,
		    median_edges,
		    q3_edges,
		    max_edges,
		    stddev_edges);

    std::cout << "[TRAVERSED EDGES SUMMARY] MEAN : " << mean_edges
	      << ", MIN : " << min_edges
	      << ", Q1 : " << q1_edges
	      << ", MEDIAN : " << median_edges
	      << ", Q3 : " << q3_edges
	      << ", MAX : " << max_edges
	      << ", STDDEV : " << stddev_edges
	      << std::endl;

#endif
  }
  

  


};

#endif
