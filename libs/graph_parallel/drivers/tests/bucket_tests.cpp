#include "../common/utils.hpp"
#include <boost/graph/agm/util/bucket.hpp>
#include <boost/tuple/tuple.hpp>

#include <limits.h>

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
    
  typedef boost::graph::agm::buckets<work_item_t, ordering> buckets_t;
  buckets_t all_buckets(_ord);

  for (int i=0; i < wis.size(); i++) {
    all_buckets.push(wis[i]);
  }

  assert(all_buckets.get_buckets_created() == expected);
  all_buckets.print_buckets();
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
  
  void operator()(int tid,
                  std::vector<work_item_t>& wis) {
    for (typename std::vector<work_item_t>::size_type i = tid ;
         i < wis.size(); i+= nthreads) {
      work_item_t& wi = wis[i];
      all_bkts.push(wi);      
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

  
  typedef boost::graph::agm::buckets<work_item_t, ordering> buckets_t;
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
  all_buckets.print_buckets();
}

template<typename ordering>
void invoke_synchronization(int argc, char* argv[],
                            int rank,
                            int world_sz,
                            ordering _ord,
                            int nthreads,
                            std::vector<work_item_t>& wis,
                            bool expected) {

  typedef boost::graph::agm::buckets<work_item_t, ordering> buckets_t;
  buckets_t all_buckets(world_sz, _ord);
  thread_executor<buckets_t> te(nthreads, all_buckets);  

  boost::scoped_array<boost::thread> threads(new boost::thread[nthreads - 1]);
  for (int i = 0; i < nthreads - 1; ++i) {
    boost::thread thr(boost::ref(te), i + 1, wis);
    threads[i].swap(thr);
  }
	  
  te(0, wis);
    
  for (int i = 0; i < (nthreads - 1); ++i)
    threads[i].join();

  if (rank == 0) { 
    std::cout << "Rank : " << rank
              << "All buckets created = " << all_buckets.get_buckets_created()
              << std::endl;
    all_buckets.print_buckets();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 1)  {
    std::cout << "Rank : " << rank
              << "All buckets created = " << all_buckets.get_buckets_created()
              << std::endl;
    all_buckets.print_buckets();
  }

  //    assert(all_buckets.get_buckets_created() == expected);

  // Barrier
  MPI_Barrier(MPI_COMM_WORLD);
  // call synchronize

  assert(all_buckets.synchronize_current_bucket() == expected);

  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "Printing buckets after barrier synchronization ..." << std::endl;

    if (rank == 0) { 
    std::cout << "[A] Rank : " << rank
              << "All buckets created = " << all_buckets.get_buckets_created()
              << std::endl;
    all_buckets.print_buckets();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 1)  {
    std::cout << "[A] Rank : " << rank
              << "All buckets created = " << all_buckets.get_buckets_created()
              << std::endl;
    all_buckets.print_buckets();
  }

}

template<typename ordering>
void test_synchronize(int argc, char* argv[],
                      int rank,
                      int world_sz,
                      ordering _ord,
                      int nthreads,
                      bool expected) {

  std::vector<work_item_t> wis;

  // create un-even buckets in both ranks
  if (rank == 0) {

    work_item_t w1(0,5);
    work_item_t w3(0,9);
    work_item_t w4(0,25);
    work_item_t w5(0,15);
    work_item_t w7(0,34);
    work_item_t w8(0,16);
    work_item_t w9(0,8);
    work_item_t w10(0,9);
    work_item_t w11(8, 45);

    wis.push_back(w1);
    wis.push_back(w3);
    wis.push_back(w4);
    wis.push_back(w5);
    wis.push_back(w7);
    wis.push_back(w8);
    wis.push_back(w9);
    wis.push_back(w10);
    wis.push_back(w11);  
  }

  if (rank == 1) {

    work_item_t w11(0,1);
    work_item_t w12(0,2);
    work_item_t w13(0,2);
    work_item_t w14(0,3);
    work_item_t w15(0,4);
    work_item_t w16(0,9);
    work_item_t w17(3,8);

    wis.push_back(w11);
    wis.push_back(w12);
    wis.push_back(w13);
    wis.push_back(w14);
    wis.push_back(w15);
    wis.push_back(w16);
    wis.push_back(w17);
  }

  invoke_synchronization(argc,
                         argv,
                         rank,
                         world_sz,
                         _ord,
                         nthreads,
                         wis,
                         expected);
}




template<typename ordering>
void test_synchronize_exchange_both(int argc, char* argv[],
                      int rank,
                      int world_sz,
                      ordering _ord,
                      int nthreads,
                      bool expected) {

  std::vector<work_item_t> wis;

  // create un-even buckets in both ranks
  if (rank == 0) {

    work_item_t w1(0,5);
    work_item_t w3(0,9);
    work_item_t w4(0,25);
    work_item_t w5(0,15);
    work_item_t w7(0,34);
    work_item_t w8(0,16);
    work_item_t w9(0,8);
    work_item_t w10(0,9);
    work_item_t w11(8, 45);

    wis.push_back(w1);
    wis.push_back(w3);
    wis.push_back(w4);
    wis.push_back(w5);
    wis.push_back(w7);
    wis.push_back(w8);
    wis.push_back(w9);
    wis.push_back(w10);
    wis.push_back(w11);  
  }

  if (rank == 1) {

    work_item_t w11(0,1);
    work_item_t w12(0,2);
    work_item_t w13(0,2);
    work_item_t w14(0,3);
    work_item_t w15(0,4);
    work_item_t w16(0,39);
    work_item_t w17(3,48);

    wis.push_back(w11);
    wis.push_back(w12);
    wis.push_back(w13);
    wis.push_back(w14);
    wis.push_back(w15);
    wis.push_back(w16);
    wis.push_back(w17);
  }

  invoke_synchronization(argc,
                         argv,
                         rank,
                         world_sz,
                         _ord,
                         nthreads,
                         wis,
                         expected);
}


template<typename ordering>
void test_synchronize_one_empty(int argc, char* argv[],
                      int rank,
                      int world_sz,
                      ordering _ord,
                      int nthreads,
                      bool expected) {

  std::vector<work_item_t> wis;

  // create un-even buckets in both ranks
  if (rank == 0) {
    // no work
  }

  if (rank == 1) {

    work_item_t w11(0,1);
    work_item_t w12(0,2);
    work_item_t w13(0,2);
    work_item_t w14(0,3);
    work_item_t w15(0,4);
    work_item_t w16(0,39);
    work_item_t w17(3,48);

    wis.push_back(w11);
    wis.push_back(w12);
    wis.push_back(w13);
    wis.push_back(w14);
    wis.push_back(w15);
    wis.push_back(w16);
    wis.push_back(w17);
  }

  invoke_synchronization(argc,
                         argv,
                         rank,
                         world_sz,
                         _ord,
                         nthreads,
                         wis,
                         expected);
}


template<typename ordering>
void test_synchronize_both_empty(int argc, char* argv[],
                      int rank,
                      int world_sz,
                      ordering _ord,
                      int nthreads,
                      bool expected) {

  std::vector<work_item_t> wis;

  // create un-even buckets in both ranks
  if (rank == 0) {
    // no work
  }

  if (rank == 1) {
    // no work
  }

  invoke_synchronization(argc,
                         argv,
                         rank,
                         world_sz,
                         _ord,
                         nthreads,
                         wis,
                         expected);
}


int main(int argc, char* argv[]) {

  MPI_Init(NULL, NULL);
  
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  _RANK = world_rank;
  
  int nthreads = 2;  
  {
    amplusplus::register_builtin_mpi_datatypes();
    amplusplus::register_mpi_datatype<work_item_t>();
    //amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item_t>());
  
    info("Testing chaotic buckets ...");
    execute<chaotic>(argc, argv, chaotic(), (uint64_t)1);

    info("Testing dijkstra buckets ...");
    execute< dijkstra<1> >(argc, argv, dijkstra<1>(), (uint64_t)10);

    info("Testing delta buckets ...");
    execute< delta_ord<1, uint64_t> >(argc, argv, delta_ord<1, uint64_t>(5), (uint64_t)6);
  
    info("Testing mt chaotic buckets ...");
    execute_mt<chaotic>(argc, argv, chaotic(), nthreads, (uint64_t)1);

    info("Testing mt dijkstra buckets ...");
    execute_mt< dijkstra<1> >(argc, argv, dijkstra<1>(), nthreads, (uint64_t)10);

    info("Testing mt delta buckets ...");
    execute_mt< delta_ord<1, uint64_t> >(argc, argv, delta_ord<1, uint64_t>(5), nthreads, (uint64_t)6);

    info("Testing synchronization ...");
    test_synchronize(argc, argv, world_rank, world_size, delta_ord<1, uint64_t>(5),
                     nthreads, false);

    info("Testing synchronization -- both inserts ...");
    test_synchronize_exchange_both(argc, argv, world_rank, world_size, delta_ord<1, uint64_t>(5),
                                   nthreads, false);

    info("Testing synchronization -- one empty ...");
    test_synchronize_one_empty(argc, argv, world_rank, world_size, delta_ord<1, uint64_t>(5),
                               nthreads, false);

    info("Testing synchronization -- both empty ...");
    test_synchronize_both_empty(argc, argv, world_rank, world_size, delta_ord<1, uint64_t>(5),
                                nthreads, true);

    
  }
  MPI_Finalize();
}
