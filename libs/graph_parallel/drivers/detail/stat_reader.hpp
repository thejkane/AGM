#ifndef PBGL2_STAT_READER
#define PBGL2_STAT_READER
#include <vector>
#include "equation_builder.hpp"

typedef std::vector<int64_t> WorkStatData_t;

template <typename graph_create_params,
	  typename instance_params,
	  typename algorithm_params,
	  typename emulator_params>
class stat_reader {

private:
  graph_create_params& ggen_params;
  instance_params& runtime_params;
  algorithm_params& algo_params;
  emulator_params& em_params;



public:
  stat_reader(graph_create_params& gcp,
	      instance_params& ip,
	      algorithm_params& ap,
	      emulator_params& ep): ggen_params(gcp),
				    runtime_params(ip),
				    algo_params(ap),
				    em_params(ep) {}

  void print(WorkStatData_t& risinguseful,
	     WorkStatData_t& risinginvalidated,
	     WorkStatData_t& risingrejected) {

    std::cout << "Diff useful : {";
    BOOST_FOREACH(uint64_t ruseful, risinguseful) {
      std::cout << ruseful << ", ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Diff invalidated(diff) : {";
    BOOST_FOREACH(uint64_t invalid, risinginvalidated) {
      std::cout << invalid << ", ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Diff rejected(diff) : {";
    BOOST_FOREACH(uint64_t reject, risingrejected) {
      std::cout << reject << ", ";
    }
    std::cout << "}" << std::endl;
  }


  void simulate(linear_equation& finvalid,
		linear_equation& freject,
		uint64_t startlevel) {
    // the parallelism is decided by the graph structure
    // for chaotic bfs. Therefore, we calculated the amount
    // of invalidated and rejected work generated in each level
    // calculates how much time-steps it takes to process the
    // generated work
    uint64_t time_steps = 0;
    for (auto i=startlevel; i < ggen_params.n; ++i) {
      uint64_t invalidatedwk = finvalid(i);
      uint64_t rejectwk = freject(i);

      std::cout << "[INFO] Level : " << i << ", invalidated : " << invalidatedwk << ", rejectwk :" << rejectwk
		<< std::endl;
      // total parallel work generated is invalidatedwk+rejectwk
      // total time to execute invalidate work
      uint64_t inv_ts = (invalidatedwk/runtime_params.threads) + 1;
      // total time to execute rejected work
      uint64_t rej_ts = (rejectwk/runtime_params.threads) + 1;

      time_steps += (inv_ts + rej_ts);
    }

    std::cout << "Total timesteps : " << time_steps << std::endl;
  }

  void process(WorkStatData_t& usefuldiff,
	       WorkStatData_t& invalidatediff,
	       WorkStatData_t& rejectdiff) {

    print(usefuldiff, invalidatediff, rejectdiff);

    uint64_t increaselevel = 0;
    // find max invalidated diff value
    uint64_t maxinvalid = 0;
    // do two times
    bool invalidfirst = true;
    BOOST_FOREACH(uint64_t invalid, invalidatediff) {
      ++increaselevel;
      if (maxinvalid < invalid) {
	maxinvalid = invalid;
	invalidfirst = true;
      } else if(maxinvalid > invalid) {
	if (invalidfirst) {
	  invalidfirst = false;
	} else {
	  // found max invalid, did for two times
	  break;
	}
      }
    }

    if (invalidfirst) {
      std::cout << "[ERROR] Level iterations are not sufficient to calculate the max invalidate derivative." << std::endl;
      return;
    }

    // find max rejected diff value
    uint64_t maxrejected = 0;
    // do two times
    bool rejectfirst = true;
    BOOST_FOREACH(uint64_t reject, rejectdiff) {
      if (maxrejected < reject) {
	maxrejected = reject;
	rejectfirst = true;
      } else if(maxrejected > reject) {
	if (rejectfirst) {
	  rejectfirst = false;
	} else {
	  // found max invalid, did for two times
	  break;
	}
      }
    }

    if (rejectfirst) {
      std::cout << "[ERROR] Level iterations are not sufficient to calculate the max rejected derivative." << std::endl;
      return;
    }
    
    std::cout << "[INFO] Maximum invalidated derivative : " << maxinvalid << ", maximum rejected derivative : " << maxrejected << std::endl;

    double minvalid = ((double)0 - (double)maxinvalid)/((double)ggen_params.n - (double)increaselevel);
    double cinvalid = (double)ggen_params.n;

    std::cout << "[INFO] The m for invalid work : " << minvalid << ", the c for invalid work : " << cinvalid << std::endl;

    linear_equation finvalid(minvalid, cinvalid);

    double mreject = ((double)0 - (double)maxrejected)/((double)ggen_params.n - (double)increaselevel);
    double creject = (double)ggen_params.n;

    std::cout << "[INFO] The m for reject work : " << mreject << ", the c for reject work : " << creject << std::endl;

    linear_equation freject(mreject, creject);

    simulate(finvalid, freject, increaselevel);
  }
};

#endif
