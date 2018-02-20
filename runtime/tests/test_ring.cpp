// Copyright 2010-2013 The Trustees of Indiana University.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#include <config.h>

#include <boost/config.hpp>
#define IS_MPI_TRANSPORT_mpi 1
#define IS_MPI_TRANSPORT_gasnet 0
#define IS_MPI_TRANSPORT_shm 0
#define IS_MPI_TRANSPORT BOOST_JOIN(IS_MPI_TRANSPORT_, TRANSPORT)
#define IS_SHM_TRANSPORT_mpi 0
#define IS_SHM_TRANSPORT_gasnet 0
#define IS_SHM_TRANSPORT_shm 1
#define IS_SHM_TRANSPORT BOOST_JOIN(IS_SHM_TRANSPORT_, TRANSPORT)

#if IS_SHM_TRANSPORT
#include <omp.h>
#endif

#include "am++/am++.hpp"
#include <boost/config.hpp>
#define TRANSPORT_HEADER <am++/BOOST_JOIN(TRANSPORT, _transport).hpp>
#include TRANSPORT_HEADER
#include "am++/basic_coalesced_message_type.hpp"
#include <stdio.h>
#include <string>
#include <numeric>
#include <unistd.h>
// #include <asm/msr.h>
// #include <valgrind/callgrind.h>
#if IS_MPI_TRANSPORT
#include <mpi.h>
#include <nbc.h>
#endif
#include <boost/assert.hpp>

inline void rdtscll(unsigned long long& x) {x = 0;}

//#define HYBRID

volatile int ctr=0;

const int msgs = 1; // 300;

typedef amplusplus::transport::rank_type rank_type;

struct ring_handler {
  typedef amplusplus::basic_coalesced_message_type<int, ring_handler> ring_msg_type;
  ring_msg_type* ring_msg;
  rank_type p,r,next;
  ring_handler() : ring_msg(NULL), p(0), r(0), next(0) {}
  ring_handler(rank_type _p, rank_type _r, ring_msg_type& ring_msg) : ring_msg(&ring_msg), p(_p), r(_r), next((_r + 1) % _p) {};
  void operator()(int msg, rank_type source) const;
};

inline void ring_handler::operator()(int msg, rank_type /*source*/) const {
  // std::cout << r << " received " << msg << " from " << source << "!\n"; 
  if(r) {(*ring_msg).message_being_built(next); (*ring_msg).send(msg, next);}

  // __sync_fetch_and_add(&ctr, 1);
  ++ctr;
};

#if IS_MPI_TRANSPORT
void run_mpi(int rank, int size) {
  // if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

  const int tests=100000;
  std::vector<double> test(tests);
  std::vector<unsigned long long> test_tsc(tests);
  size_t datasize = 32 * msgs;
  char* data = (char*)alloca(datasize);
  BOOST_ASSERT (data);
  for(int i=0; i<tests; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);

    test[i] = -MPI_Wtime();
    rdtscll(test_tsc[i]);
    test_tsc[i] = -test_tsc[i];

    if (rank == 0) {
      MPI_Send(data, datasize, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
      MPI_Recv(data, datasize, MPI_BYTE, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(data, datasize, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(data, datasize, MPI_BYTE, (rank + 1) % size, 0, MPI_COMM_WORLD);
    }

    test[i] += MPI_Wtime();
    test[i]*=1e6;
    unsigned long long tsc2;
    rdtscll(tsc2);
    test_tsc[i] += tsc2;
  }
  double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
  double min = *min_element(test.begin(), test.end());
  double max = *max_element(test.begin(), test.end());
  std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
  nth_element(test.begin(), nthblock, test.end());
  double med = *nthblock;
  unsigned long long avg_tsc = std::accumulate(test_tsc.begin(), test_tsc.end(), (unsigned long long)0)/(unsigned long long)test_tsc.size();
  unsigned long long min_tsc = *min_element(test_tsc.begin(), test_tsc.end());
  unsigned long long max_tsc = *max_element(test_tsc.begin(), test_tsc.end());
  std::vector<unsigned long long>::iterator nthblock_tsc = test_tsc.begin()+test_tsc.size()/2;
  nth_element(test_tsc.begin(), nthblock_tsc, test_tsc.end());
  unsigned long long med_tsc = *nthblock_tsc;
  if(!rank) std::cout << "MPI (in us) min: " << min << " max: " << max <<" avg: "<<avg<<" med: "<<med<<"\n";
  if(!rank) std::cout << "MPI (in ticks) min: " << min_tsc << " max: " << max_tsc <<" avg: "<<avg_tsc<<" med: "<<med_tsc<<"\n";

}

void run_mpi2os(int rank, int size) {
  if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

  char remdata=255;
  char data = (rank+1)%size;

  MPI_Win win;
  MPI_Win_create(&remdata, sizeof(remdata), 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  const int tests=1000;
  std::vector<double> test(tests);

  for(int i=0; i<tests; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);

    test[i] = -MPI_Wtime();

    if (rank == 0) {

      //printf("[%i] %i %i\n", rank, (int)remdata, (int)data);
      
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 1, 0, win);
      MPI_Put(&data, 1, MPI_CHAR, 1, 0, 1, MPI_CHAR, win);
      MPI_Win_unlock(1, win);
      
      //printf("[%i] %i %i\n", rank, (int)remdata, (int)data);
      
      while(remdata != rank) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, win);
        MPI_Win_unlock(rank, win);
      };

    } else {

      //printf("[%i] %i %i\n", rank, (int)remdata, (int)data);

      while(remdata != rank) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, win);
        MPI_Win_unlock(rank, win);
      };

      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, (rank+1)%size, 0, win);
      MPI_Put(&data, 1, MPI_CHAR, (rank+1)%size, 0, 1, MPI_CHAR, win);
      MPI_Win_unlock((rank+1)%size, win);

      //printf("[%i] %i %i\n", rank, (int)remdata, (int)data);

    }

    test[i] += MPI_Wtime();
    test[i]*=1e6;
  }
  double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
  double min = *min_element(test.begin(), test.end());
  double max = *max_element(test.begin(), test.end());
  std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
  nth_element(test.begin(), nthblock, test.end());
  double med = *nthblock;
  if(!rank) std::cout << "MPI2OS (in us) min: " << min << " max: " << max <<" avg: "<<avg<<" med: "<<med<<"\n";

  MPI_Win_free(&win);
}


void run_mpi_alloc(int rank, int size) {
  if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

  const int tests=100000;
  std::vector<double> test(tests);
  std::vector<unsigned long long> test_tsc(tests);
  size_t datasize = 32 * msgs;
  char* data;
  for(int i=0; i<tests; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);

    test[i] = -MPI_Wtime();
    rdtscll(test_tsc[i]);
    test_tsc[i] = -test_tsc[i];

    MPI_Alloc_mem(datasize, MPI_INFO_NULL, (void*)&data);
    if (rank == 0) {
      MPI_Send(data, datasize, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
      MPI_Recv(data, datasize, MPI_BYTE, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(data, datasize, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(data, datasize, MPI_BYTE, (rank + 1) % size, 0, MPI_COMM_WORLD);
    }

    MPI_Free_mem(data);
    test[i] += MPI_Wtime();
    test[i]*=1e6;
    unsigned long long tsc2;
    rdtscll(tsc2);
    test_tsc[i] += tsc2;
  }
  double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
  double min = *min_element(test.begin(), test.end());
  double max = *max_element(test.begin(), test.end());
  std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
  nth_element(test.begin(), nthblock, test.end());
  double med = *nthblock;
  unsigned long long avg_tsc = std::accumulate(test_tsc.begin(), test_tsc.end(), (unsigned long long)0)/(unsigned long long)test_tsc.size();
  unsigned long long min_tsc = *min_element(test_tsc.begin(), test_tsc.end());
  unsigned long long max_tsc = *max_element(test_tsc.begin(), test_tsc.end());
  std::vector<unsigned long long>::iterator nthblock_tsc = test_tsc.begin()+test_tsc.size()/2;
  nth_element(test_tsc.begin(), nthblock_tsc, test_tsc.end());
  unsigned long long med_tsc = *nthblock_tsc;
  if(!rank) std::cout << "MPI alloc (in us) min: " << min << " max: " << max <<" avg: "<<avg<<" med: "<<med<<"\n";
  if(!rank) std::cout << "MPI alloc (in ticks) min: " << min_tsc << " max: " << max_tsc <<" avg: "<<avg_tsc<<" med: "<<med_tsc<<"\n";

}

void run_mpi_test(int rank, int size) {
  if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

  const int tests=100000;
  std::vector<double> test(tests);
  std::vector<unsigned long long> test_tsc(tests);
  size_t datasize = 32 * msgs;
  char* data = (char*)alloca(datasize);
  BOOST_ASSERT (data);
  for(int i=0; i<tests; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);

    test[i] = -MPI_Wtime();
    rdtscll(test_tsc[i]);
    test_tsc[i] = -test_tsc[i];

    MPI_Status stat;
    if (rank == 0) {
      MPI_Send(data, datasize, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
      MPI_Request req;
      MPI_Irecv(data, datasize, MPI_BYTE, size - 1, 0, MPI_COMM_WORLD, &req);
      while (true) {
        int flag = 0;
        MPI_Test(&req, &flag, &stat);
        if (flag) {
          break;
        }
      }
    } else {
      MPI_Request req;
      MPI_Irecv(data, datasize, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, &req);
      while (true) {
        int flag = 0;
        MPI_Test(&req, &flag, &stat);
        if (flag) {
          MPI_Send(data, datasize, MPI_BYTE, (rank + 1) % size, 0, MPI_COMM_WORLD);
          break;
        }
      }
    }

    test[i] += MPI_Wtime();
    test[i]*=1e6;
    unsigned long long tsc2;
    rdtscll(tsc2);
    test_tsc[i] += tsc2;
  }
  double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
  double min = *min_element(test.begin(), test.end());
  double max = *max_element(test.begin(), test.end());
  std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
  nth_element(test.begin(), nthblock, test.end());
  double med = *nthblock;
  unsigned long long avg_tsc = std::accumulate(test_tsc.begin(), test_tsc.end(), (unsigned long long)0)/(unsigned long long)test_tsc.size();
  unsigned long long min_tsc = *min_element(test_tsc.begin(), test_tsc.end());
  unsigned long long max_tsc = *max_element(test_tsc.begin(), test_tsc.end());
  std::vector<unsigned long long>::iterator nthblock_tsc = test_tsc.begin()+test_tsc.size()/2;
  nth_element(test_tsc.begin(), nthblock_tsc, test_tsc.end());
  unsigned long long med_tsc = *nthblock_tsc;
  if(!rank) std::cout << "MPI Test (in us) min: " << min << " max: " << max <<" avg: "<<avg<<" med: "<<med<<"\n";
  if(!rank) std::cout << "MPI Test (in ticks) min: " << min_tsc << " max: " << max_tsc <<" avg: "<<avg_tsc<<" med: "<<med_tsc<<"\n";

}

void run_mpi_iprobe(int rank, int size) {
  if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

  const int tests=100000;
  std::vector<double> test(tests);
  std::vector<unsigned long long> test_tsc(tests);
  size_t datasize = 32 * msgs;
  char* data = (char*)alloca(datasize);
  BOOST_ASSERT (data);
  for(int i=0; i<tests; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);

    test[i] = -MPI_Wtime();
    rdtscll(test_tsc[i]);
    test_tsc[i] = -test_tsc[i];

    MPI_Status stat;
    if (rank == 0) {
      MPI_Send(data, datasize, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
      while (true) {
        int flag = 0;
        MPI_Iprobe(size - 1, 0, MPI_COMM_WORLD, &flag, &stat);
        if (flag) {
          MPI_Recv(data, datasize, MPI_BYTE, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          break;
        }
      }
    } else {
      while (true) {
        int flag = 0;
        MPI_Iprobe(rank - 1, 0, MPI_COMM_WORLD, &flag, &stat);
        if (flag) {
          MPI_Recv(data, datasize, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Send(data, datasize, MPI_BYTE, (rank + 1) % size, 0, MPI_COMM_WORLD);
          break;
        }
      }
    }

    test[i] += MPI_Wtime();
    test[i]*=1e6;
    unsigned long long tsc2;
    rdtscll(tsc2);
    test_tsc[i] += tsc2;
  }
  double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
  double min = *min_element(test.begin(), test.end());
  double max = *max_element(test.begin(), test.end());
  std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
  nth_element(test.begin(), nthblock, test.end());
  double med = *nthblock;
  unsigned long long avg_tsc = std::accumulate(test_tsc.begin(), test_tsc.end(), (unsigned long long)0)/(unsigned long long)test_tsc.size();
  unsigned long long min_tsc = *min_element(test_tsc.begin(), test_tsc.end());
  unsigned long long max_tsc = *max_element(test_tsc.begin(), test_tsc.end());
  std::vector<unsigned long long>::iterator nthblock_tsc = test_tsc.begin()+test_tsc.size()/2;
  nth_element(test_tsc.begin(), nthblock_tsc, test_tsc.end());
  unsigned long long med_tsc = *nthblock_tsc;
  if(!rank) std::cout << "MPI Iprobe (in us) min: " << min << " max: " << max <<" avg: "<<avg<<" med: "<<med<<"\n";
  if(!rank) std::cout << "MPI Iprobe (in ticks) min: " << min_tsc << " max: " << max_tsc <<" avg: "<<avg_tsc<<" med: "<<med_tsc<<"\n";

}
#endif

#if 0
void run_am(int rank, int size) {
  { // am needs to be destructed before MPI_Finalize()
    const int tests=1000000;
    std::vector<double> test(tests);
    std::vector<unsigned long long> test_tsc(tests);
    if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

    // CALLGRIND_START_INSTRUMENTATION
    am_type am;
    boost::message_sender<am_type, ring_handler(int)> send_ring(am);
    send_ring.set_handler(ring_handler(size, rank, send_ring));

    for(int i=0; i<tests; ++i) {
      // MPI_Barrier(MPI_COMM_WORLD);
      am.barrier();

      test[i] = -MPI_Wtime();
      rdtscll(test_tsc[i]);
      test_tsc[i] = -test_tsc[i];

      if (!rank) for(int i=0; i<msgs; ++i) { send_ring[1](i); }
      am.flush();
      // fprintf(stderr, "All messages sent\n");
      
      //am.synchronize();
      //am.quiesce();
      while(msgs != ctr) {am.progress();}
      // am.synchronize();
      am.flush();
      // fprintf(stderr, "Done looping\n");
      ctr = 0;
      // while(msgs != __sync_val_compare_and_swap(&ctr, msgs, 0));

      test[i] += MPI_Wtime();
      test[i]*=1e6;
      unsigned long long tsc2;
      rdtscll(tsc2);
      test_tsc[i] += tsc2;
    }
    // CALLGRIND_STOP_INSTRUMENTATION
    double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
    double min = *min_element(test.begin(), test.end());
    double max = *max_element(test.begin(), test.end());
    std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
    nth_element(test.begin(), nthblock, test.end());
    double med = *nthblock;
    unsigned long long avg_tsc = std::accumulate(test_tsc.begin(), test_tsc.end(), (unsigned long long)0)/(unsigned long long)test_tsc.size();
    unsigned long long min_tsc = *min_element(test_tsc.begin(), test_tsc.end());
    unsigned long long max_tsc = *max_element(test_tsc.begin(), test_tsc.end());
    std::vector<unsigned long long>::iterator nthblock_tsc = test_tsc.begin()+test_tsc.size()/2;
    nth_element(test_tsc.begin(), nthblock_tsc, test_tsc.end());
    unsigned long long med_tsc = *nthblock_tsc;
    if(!rank) std::cout << "AM++ (in us) min: " << min << " max: " << max <<" avg: "<<avg<<" med: "<<med<<"\n";
    if(!rank) std::cout << "AM++ (in ticks) min: " << min_tsc << " max: " << max_tsc <<" avg: "<<avg_tsc<<" med: "<<med_tsc<<"\n";

  }
}
#endif

void run_am(amplusplus::environment& env, rank_type rank, rank_type size) {
  amplusplus::transport transport = env.create_transport();
  amplusplus::basic_coalesced_message_type<int, ring_handler> ring_msg(amplusplus::basic_coalesced_message_type_gen(1 << 14), transport);
  ring_msg.set_handler(ring_handler(size, rank, ring_msg));
  { // am needs to be destructed before MPI_Finalize()
    const int tests=1000;
    std::vector<double> test(tests);
    std::vector<unsigned long long> test_tsc(tests);
    // if (size < 2) {fprintf(stderr, "Need at least two processes\n"); abort();}

    // CALLGRIND_START_INSTRUMENTATION
    for(int i=0; i<tests; ++i) {
      // if (!rank && i % 1 == 0) std::cerr << "Test " << i + 1 << " of " << tests << std::endl;
      {amplusplus::scoped_epoch e(transport);}

      test[i] = -amplusplus::get_time();
      rdtscll(test_tsc[i]);
      test_tsc[i] = -test_tsc[i];

      {
        amplusplus::scoped_epoch epoch(transport);
        if (!rank) for(int i=0; i<msgs; ++i) {ring_msg.message_being_built(1 % size); ring_msg.send(i, 1 % size);}
      }
      ctr = 0;

      test[i] += amplusplus::get_time();
      test[i]*=1e6;
      unsigned long long tsc2;
      rdtscll(tsc2);
      test_tsc[i] += tsc2;

      {amplusplus::scoped_epoch e(transport);}
    }
    // CALLGRIND_STOP_INSTRUMENTATION
    double avg = std::accumulate(test.begin(), test.end(), (double)0)/(double)test.size();
    double mn = *min_element(test.begin(), test.end());
    double mx = *max_element(test.begin(), test.end());
    std::vector<double>::iterator nthblock = test.begin()+test.size()/2;
    nth_element(test.begin(), nthblock, test.end());
    double med = *nthblock;
    unsigned long long avg_tsc = std::accumulate(test_tsc.begin(), test_tsc.end(), (unsigned long long)0)/(unsigned long long)test_tsc.size();
    unsigned long long min_tsc = *min_element(test_tsc.begin(), test_tsc.end());
    unsigned long long max_tsc = *max_element(test_tsc.begin(), test_tsc.end());
    std::vector<unsigned long long>::iterator nthblock_tsc = test_tsc.begin()+test_tsc.size()/2;
    nth_element(test_tsc.begin(), nthblock_tsc, test_tsc.end());
    unsigned long long med_tsc = *nthblock_tsc;
    if(!rank) std::cout << "AM++ (in us) min: " << mn << " max: " << mx <<" avg: "<<avg<<" med: "<<med<<"\n";
    if(!rank) std::cout << "AM++ (in ticks) min: " << min_tsc << " max: " << max_tsc <<" avg: "<<avg_tsc<<" med: "<<med_tsc<<"\n";

  }
}

void do_one_thread(amplusplus::environment& env) {
  amplusplus::rank_type rank, size;

  {
    amplusplus::transport trans = env.create_transport();
    rank = trans.rank();
    size = trans.size();
  }

#if IS_MPI_TRANSPORT
  NBC_Request nbcreq;
  NBC_Ibarrier(MPI_COMM_WORLD, &nbcreq);
  NBC_Wait(&nbcreq, MPI_STATUS_IGNORE);
#endif

  // run_mpi(rank, size);
  // run_mpi_alloc(rank, size);
  // run_mpi_test(rank, size);
  // run_mpi_iprobe(rank, size);

  // run_mpi2os(rank, size);

  // fprintf(stderr, "Start profiling now %d\n", getpid()); fflush(stderr);
  // sleep(10);
  fprintf(stderr, "Running test\n"); fflush(stderr);
  run_am(env, rank, size);
  fprintf(stderr, "Done testing\n"); fflush(stderr);
}

int main(int argc, char** argv) {
  (void)argc; (void)argv;
#if IS_MPI_TRANSPORT
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  do_one_thread(env);
#elif IS_SHM_TRANSPORT
  boost::scoped_ptr<amplusplus::shm_environment_common> common;

#pragma omp parallel
  {
#pragma omp single
    {
      common.reset(new amplusplus::shm_environment_common(omp_get_num_threads()));
    }
    amplusplus::environment env = amplusplus::shm_environment(*common, omp_get_thread_num());
    do_one_thread(env);
  }
#endif
  return 0;
}
