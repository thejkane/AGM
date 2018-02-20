#include <config.h>
#include <boost/config.hpp>
#include "am++/am++.hpp"
#include <iostream>
#include <pthread.h>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/bind.hpp>
#define TRANSPORT_HEADER <am++/BOOST_JOIN(TRANSPORT, _transport).hpp>
#include TRANSPORT_HEADER

#include <time.h>
#include <sys/time.h>

#define DISABLE_SELF_SEND_CHECK

typedef double time_type;
// copies directly from europar_tests
inline time_type get_time()
{
  return MPI_Wtime();
#if 0
  timeval tp;
  gettimeofday(&tp, 0);
  return tp.tv_sec + tp.tv_usec / 1000000.0;
#endif
}

time_type non_priority_time = 0;
time_type priority_time = 0;
time_type time_at_start = 0;

struct message_priority_handler {
  typedef amplusplus::message_type<int> tm_type;
  tm_type& tm;
  amplusplus::transport& trans;
  message_priority_handler(tm_type& tm, amplusplus::transport& t): tm(tm), trans(t) { std::cerr << "Constructor" << std::endl; }
  void operator()(int source, const int* buf, int count) const {
//  	std::cout << "My rank is " << trans.rank() << " Messag is from rank " << source << std::endl;
	for (int i = 0; i < count; ++i) {
	  // 	  std::cout << "The message content - " << buf[i] << std::endl;
	  if (buf[i] % 2 == 0) {
		non_priority_time = get_time(); // even - non priority message
	  } else {
		priority_time = get_time(); // odd - priority message
	  }
	}
  }
};

struct empty_deleter {
  typedef void result_type;
  void operator()() const {}
  template <typename T> void operator()(const T&) const {}
};

#define MSG_SIZE 10000
struct msg_send_args {
  amplusplus::transport& trans;
  amplusplus::message_type<int>& tm;
  amplusplus::message_type<int>& priority_tm;
  bool coalesced;
  int start;
  boost::barrier& barier;

  msg_send_args(amplusplus::transport& t, amplusplus::message_type<int>& mt, amplusplus::message_type<int>& priority_m, bool c, int st,
				boost::barrier& cur_barier):
	trans(t), tm(mt), priority_tm(priority_m), coalesced(c), start(st), barier(cur_barier) {}
};


void *send_normal_messages(void *arguments) {
 
  struct msg_send_args *args = (struct msg_send_args *)arguments;
  amplusplus::transport& trans = args->trans;
  amplusplus::message_type<int> tm = args->tm;
  amplusplus::message_type<int> priority_tm = args->priority_tm;
  
  if (!(args->coalesced)) {

	// wait till both threads reach here
	args->barier.wait();
	
	// message sending must be within an epoch
	amplusplus::scoped_epoch epoch(trans);
	{
	  if (trans.rank() == 0) {
		// boost::shared_ptr<int> msg = boost::make_shared<int>(0);
		for(int j=0; j<MSG_SIZE; ++j) {
		  tm.message_being_built(1);
		  priority_tm.message_being_built(1);
		  
		  tm.send(&j, 1, 1, empty_deleter()); // even - non priority
		  ++j;
		  priority_tm.send(&j, 1, 1, empty_deleter()); //odd - priority
		}
	  }
	}
  } else {
	// TODO
  }

  return NULL;
}


void run_test(amplusplus::environment& env) {

  std::cout << "Inside run test" << std::endl;
  amplusplus::transport trans = env.create_transport();
  std::cout << "Transport size " << trans.size() << std::endl;
  
  BOOST_ASSERT (trans.size() == 2);

  // sending normal messages
  amplusplus::message_type<int> tm = trans.create_message_type<int>();
  tm.set_max_count(1);
  tm.set_handler(message_priority_handler(tm, trans));

  // sending priority messages
  amplusplus::message_type<int> ptm = trans.create_message_type<int>(1);
  ptm.set_max_count(1);
  ptm.set_handler(message_priority_handler(ptm, trans));
  
  bool coalesced = false;

  // create the barier for 2 threads
  boost::barrier bar(2);
  
  struct msg_send_args args(trans,tm, ptm, coalesced, 0, bar);
//   struct msg_send_args pargs(trans,ptm, coalesced, MSG_SIZE);

  if (trans.rank() == 0) {
	// set number of threads to 2
 	trans.set_nthreads(2);

//	{ amplusplus::scoped_epoch epoch(trans); }
	
	// create threads
	pthread_t normal_thread, priority_thread;
	int ret = 0;

	// create thread to send priority messages
	ret = pthread_create(&priority_thread, NULL, send_normal_messages, (void*)&args);
	if(ret) {
	  std::cerr << "ERROR - Error creating the normal thread. Error code - " << ret << std::endl;
	  return;
	}
	
	// create thread to send non priority messages
	ret = pthread_create(&normal_thread, NULL, send_normal_messages, (void*)&args);
	if(ret) {
	  std::cerr << "ERROR - Error creating the normal thread. Error code - " << ret << std::endl;
	  return;
	}

	
	// wait till threads finishes
	pthread_join(normal_thread, NULL);
	pthread_join(priority_thread, NULL);

	trans.set_nthreads(1);
	std::cout << "Finish sending messages ..." << std::endl;
	
  } else {
	// do nothing, but need epoch for termination detection
	amplusplus::scoped_epoch epoch(trans);

	// register handler for the message types
  }

  { amplusplus::scoped_epoch epoch(trans); }

  if (trans.rank() == 1) {
	non_priority_time = non_priority_time - time_at_start;
	priority_time = priority_time - time_at_start;
	std::cout << "The priority time - " << priority_time << " non-priority time - " << non_priority_time << std::endl;
	// for rank 1 - we should be receiving priority messages before non priority messages
	std::cout << "The difference - " << non_priority_time - priority_time <<  " the ratio - " << non_priority_time / priority_time <<  std::endl;
	
	BOOST_ASSERT(priority_time < non_priority_time * 0.99);
  }
  
}


int main(int argc, char** argv) {
  std::cout << "Main starting ..." << std::endl;
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv, false, 1, 10000);
  time_at_start = get_time();
  run_test(env);


}
