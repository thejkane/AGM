#include <iostream>
#include <time.h> 
#include <chrono>
#include <utility>
#include <am++/am++.hpp>
#include <boost/lexical_cast.hpp>
#include <mpi.h>
#include <am++/mpi_transport.hpp>

typedef std::pair<uint64_t, uint64_t> work_item_t;

int _RANK;

class params {
public:
  int recev_depth;
  int polls;
  int flow_control;
  int iterations;
  bool formated;

  params() : recev_depth(1),
	     polls(1),
	     flow_control(10),
	     iterations(1),
	     formated(false){}

  bool parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--recev-depth") {
	recev_depth = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--polls") {
	polls = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--flow-control") {
	flow_control = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--iterations") {
	iterations = boost::lexical_cast<int>( argv[i+1]);
      }

      if (arg == "--format") {
	formated = true;
      }
    }

    return true;
  }

  void print() {
    if (_RANK == 0) {
      std::cout << "Receive depth : " << recev_depth << std::endl;
      std::cout << "Polls : " << polls << std::endl;
      std::cout << "Flow Control : " << flow_control << std::endl;
      std::cout << "Iterations : " << iterations << std::endl;
    }
  }
};


struct empty_deleter {
private:
  work_item_t* pw;
public:
  typedef void result_type;
  empty_deleter(work_item_t* p) : pw(p) {}
  void operator()() const {
    delete pw;
  }
};


struct reduce_int_deleter {
private:
  uint64_t* pt;
public:
  typedef void result_type;
  reduce_int_deleter(uint64_t* p) : pt(p) {}
  void operator()() const {
    delete pt;
  }
};


struct reduce_handler {
  std::vector<uint64_t>& iteration_times;
  amplusplus::transport& trans;

  reduce_handler(amplusplus::transport& t,
		 std::vector<uint64_t>& itetimes) : trans(t),
						    iteration_times(itetimes){}
  void operator()(int source, uint64_t* v, int count) const {
    assert(trans.rank() < source);
    iteration_times[source] += *v;
  }
};

struct overhead_handler {

#ifdef CYCLE_TYPE
  typedef clock_t time_point_t;
#else
  typedef std::chrono::high_resolution_clock::time_point time_point_t;
#endif

  typedef amplusplus::message_type<work_item_t> tm_type;
  tm_type& tm;
  amplusplus::transport& trans;
  std::vector<time_point_t>& rank_times;
  amplusplus::transport::rank_type originate_source;

  overhead_handler(tm_type& tm,
		   amplusplus::transport& t,
		   std::vector<time_point_t>& rt,
		   amplusplus::transport::rank_type s): tm(tm),
							trans(t),
							rank_times(rt),
							originate_source(s){}

  void operator()(int source, const work_item_t* pw, int count) const {
    if (source == originate_source) { // rank i getting a message from rank 0, send it back
      if (trans.rank() != originate_source) {
	work_item_t* pw = new work_item_t(1, 1);
	tm.message_being_built(source);
	amplusplus::transport::rank_type s = originate_source;
#ifdef PRINT_DEBUG      
	std::cout << "Rank : " << trans.rank() << "received from 0" << std::endl;
#endif
	tm.send(pw, (size_t)1, s, empty_deleter(pw));      
      } else {
	// record time for rank 0
#ifdef CYCLE_TYPE
        rank_times[originate_source] = clock();
#else
	rank_times[originate_source] = std::chrono::high_resolution_clock::now();
#endif
      }
    } else { // rank 0 getting the reply back from other ranks
      // need to time this value
      if (trans.rank() == originate_source) {
#ifdef PRINT_DEBUG      
	std::cout << "Rank : 0 received from : " << source << std::endl; 
#endif
	// record times for other ranks
#ifdef CYCLE_TYPE
	rank_times[source] = clock();
#else
	rank_times[source] = std::chrono::high_resolution_clock::now();
#endif
      }
    }
  }
};

class overhead_calculator {

#ifdef CYCLE_TYPE
  typedef clock_t time_point_t;
#else
  typedef std::chrono::high_resolution_clock::time_point time_point_t;
#endif

private:
  time_point_t source_start;
  std::vector<time_point_t> rank_times;
  amplusplus::transport& trans;
  amplusplus::transport::rank_type source;
  int iterations;
  std::vector<uint64_t> iteration_times;

public:
  overhead_calculator(amplusplus::transport& t,
		      amplusplus::transport::rank_type s,
		      int iter) : trans(t),
				  source(s),
				  iterations(iter){
    iteration_times.resize(trans.size(), 0);
  }

  void calculate_overhead() {
    rank_times.resize(trans.size());

    typedef amplusplus::message_type<work_item_t> tm_type;
    tm_type tm = trans.create_message_type<work_item_t>();

    tm.set_max_count(1);
    tm.set_handler(overhead_handler(tm, trans, rank_times, source));

    for (int j=0; j < iterations; ++j) {
      {
	amplusplus::scoped_epoch epoch(trans);
	if (trans.rank() == source) {
	  amplusplus::transport::rank_type i;
	  for (i=0; i < trans.size(); ++i) {
#ifdef PRINT_DEBUG      
	    std::cout << "Rank 0 sending to : " << i << std::endl;
#endif
	    work_item_t* pw = new work_item_t(1, 1);
	    tm.message_being_built(i);
#ifdef CYCLE_TYPE
            source_start = clock();
#else
	    source_start = std::chrono::high_resolution_clock::now();
#endif
	    tm.send(pw, 1, i, empty_deleter(pw));
	  }
	}
      }
      
      if (_RANK == source) {
	for(amplusplus::transport::rank_type k=0; k < trans.size(); ++k) {
#ifdef CYCLE_TYPE
	  auto t = rank_times[k]-source_start;
#else
	  auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(rank_times[k]-source_start).count();
#endif
	  //	  std::cout << source << "x --> " << k << " : " << (t/2) << " (nsec)" << std::endl; 
	  iteration_times[k] += (t/2);
	}
      }

      { amplusplus::scoped_epoch epoch(trans); }

    }
  }

  void reduce_values() {
    typedef amplusplus::message_type<uint64_t> tm_type;
    tm_type tm = trans.create_message_type<uint64_t>();

    tm.set_max_count(1);
    tm.set_handler(reduce_handler(trans, iteration_times));

    {
      amplusplus::scoped_epoch epoch(trans);
      for(int i=(trans.rank()-1); i >= 0; --i) {	
	tm.message_being_built(i);
	uint64_t* pt = new uint64_t;
	*pt = iteration_times[i];
	tm.send(pt, 1, i, reduce_int_deleter(pt));
      }
    }
  }


  void print(bool formated=false) {
    if (_RANK == source) {
      for (int i=0; i < trans.size(); ++i) {
	if (formated)
	  std::cout << source << "\t" << i << "\t" << (iteration_times[i]/iterations) << std::endl;
	else
	  std::cout << source << " --> " << i << " : " << (iteration_times[i]/iterations) << " (nsec)" << std::endl; 
      }
    }
  }
};


int main(int argc, char* argv[]) {

  std::cout << "Starting latency callculator ..." << std::endl;

  params p;
  p.parse(argc, argv);

  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, true /*need threading*/,
							    p.recev_depth,
							    p.polls,
							    p.flow_control);
  amplusplus::transport trans = env.create_transport();
  MPI_Comm_rank(MPI_COMM_WORLD, &_RANK);

  p.print();

  amplusplus::register_mpi_datatype<work_item_t>();

  for (int i=0; i < trans.size(); ++i) {
    overhead_calculator oc(trans, i, p.iterations);
    oc.calculate_overhead();
    //    oc.reduce_values();
    { amplusplus::scoped_epoch epoch(trans); }
    oc.print(p.formated);
  }

  return 0;
}
