// Copyright (C) 2017 The Trustees of Indiana University.
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Triangle Counting Algortihms================//
// Driver for Triangle Counting (TC) family of algorithms.
//===========================================================//

#include <iostream>
#include "common/synthetic_generator.hpp"
#include "common/parser.hpp"
#include "common/executor.hpp"
#include "common/vertex_permutations.hpp"
#include "common/spinlock.hpp"
#include <boost/graph/distributed/triangle_counting.hpp>
#include <boost/graph/distributed/triangle_counting_blocked.hpp>
#include <boost/graph/distributed/triangle_counting_sucsuc.hpp>
#include <boost/graph/distributed/triangle_counting_sps.hpp>
#include <boost/graph/distributed/triangle_counting_level.hpp>

enum tc_algorithm {
  dag_tc_striped,
  dag_tc_blocked,
  dag_tc_sps,
  dag_tc_level,
  dag_tc_predpred,
  dag_tc_sucsuc
};


class tc_instance_params {

public:
  tc_algorithm algorithm;
  id_distribution_t id_distribution;
  bool verify;
  uint64_t triangles;
  uint64_t block_size = 1000;
  uint64_t suc_block_size = 2;
  bool degreeasc;

  tc_instance_params(tc_algorithm& alg,
		     id_distribution_t& idd,
		     bool v,
		     uint64_t& t,
		     uint64_t bz,
		     uint64_t sucbz,
		     bool degord):
    algorithm(alg), id_distribution(idd), verify(v), triangles(t), block_size(bz), 
    suc_block_size(sucbz),
    degreeasc(degord){}

  void print() {

    std::cout << "Algorithm : " << get_algorithm() << std::endl;

    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == degree) {
      if (degreeasc)
	std::cout << "id distribution : degree-asc" << std::endl;
      else
	std::cout << "id distribution : degree-dsc" << std::endl;
    }

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "verify : " << verify << std::endl;
    std::cout << "block size : " << block_size << std::endl;
    std::cout << "successor block size : " << suc_block_size << std::endl;
    std::cout << "expected triangles : " << triangles << std::endl;
  }

  std::string get_algorithm() {
    std::string ordered = "-Ordered";
#ifdef TC_NO_ORDERING
    ordered = "";
#endif
    if (algorithm == dag_tc_striped) 
      return ("TC Striped"+ordered);
    else if(algorithm == dag_tc_blocked)
      return ("TC Blocked(PSP)"+ordered);
    else if(algorithm == dag_tc_sps)
      return ("TC SPS"+ordered);
    else if(algorithm == dag_tc_level)
      return ("TC Level"+ordered);
    else if(algorithm == dag_tc_predpred)
      return ("TC PredPred"+ordered);
    else if(algorithm == dag_tc_sucsuc) {
#ifdef TC_ALGO_SUCSUC
      return ("TC SucSuc"+ordered);
#else
      return ("TC PredPred"+ordered);
#endif
    } else
      return "Invalid Algorithm!";
  }
};


class tc_params {

private:
  std::vector<tc_instance_params> instance_params;
  std::vector<tc_algorithm> algorithms;
  id_distribution_t id_distribution = horizontal; // default
  bool verify = false;
  uint64_t triangles;
  std::vector<uint64_t> block_size;
  std::vector<uint64_t> suc_block_sizes;
  bool degreeasc = true;

public:
  bool parse(int argc, char* argv[]){
    for (int i = 1; i < argc; ++i) {
      tc_algorithm algorithm;
      std::string arg = argv[i];
      if (arg == "--run_tc_striped") {
	algorithm = dag_tc_striped;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run_tc_blocked") {
	algorithm = dag_tc_blocked;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run_tc_sps") {
	algorithm = dag_tc_sps;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run_tc_sucsuc") {
	algorithm = dag_tc_sucsuc;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run_tc_level") {
	algorithm = dag_tc_level;
	algorithms.push_back(algorithm);
      }

      if (arg == "--pred-block-size") {
	block_size = extract_params<size_t>( argv[i+1] );
      }

      if (arg == "--suc-block-size") {
	suc_block_sizes = extract_params<size_t>( argv[i+1] );
      }

      if (arg == "--expected-triangles")
	triangles = boost::lexical_cast<uint64_t>( argv[i+1]);

      if (arg == "--id-distribution") {
	if (strcmp(argv[i+1],"vertical") == 0)
	  id_distribution = vertical;
	else if (strcmp(argv[i+1],"horizontal") == 0)
	  id_distribution = horizontal;
	else if (strcmp(argv[i+1],"degree-asc") == 0) {
	  id_distribution = degree;
	  degreeasc = true;
	} else if (strcmp(argv[i+1],"degree-dsc") == 0) {
	  id_distribution = degree;
	  degreeasc = false;
	} else {
	  std::cout << "Invalid id distribution type. Available types are vertical and horizontal" 
		    << std::endl;
	  return false;
	}
      }

      if (arg == "--verify") {
	verify = true;
      }
    }

    if (suc_block_sizes.empty())
      suc_block_sizes.push_back(2);

    if (block_size.empty())
      block_size.push_back(1000);
  }

  void print() {
    BOOST_FOREACH(tc_algorithm alg, algorithms) {
      std::cout << "Algorithm -- ";
      if (alg == dag_tc_striped) 
	std::cout << "TC Synchronous";
      else if(alg == dag_tc_blocked)
	std::cout << "TC PSP";
      else if(alg == dag_tc_sps)
	std::cout << "TC SPS";
      else if(alg == dag_tc_predpred)
	std::cout << "TC Predecessor-Predecessor";
      else if(alg == dag_tc_sucsuc) {
#ifdef TC_ALGO_SUCSUC
	std::cout << "TC Successor-Successor";
#else
	std::cout << "TC Predecessor-Predecessor";
#endif
      } else
	std::cerr << "Invalid !" << std::endl;

      std::cout << std::endl;
    }

    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == degree) {
      if (degreeasc)
	std::cout << "id distribution : degree-asc" << std::endl;
      else
	std::cout << "id distribution : degree-dsc" << std::endl;
    }
    
    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "verify : " << verify << std::endl;
    std::cout << "expected triangles : " << triangles << std::endl;
  }

  std::vector<tc_instance_params>&
  get_instance_params() {
    if (instance_params.empty()) {
      BOOST_FOREACH(uint64_t bsz, block_size) { 
	BOOST_FOREACH(uint64_t sucbsz, suc_block_sizes) { 
	  BOOST_FOREACH(tc_algorithm alg, algorithms) {
	    instance_params.push_back(tc_instance_params(alg, id_distribution, 
							 verify, triangles,
							 bsz, sucbsz, degreeasc));
	  }
	}
      }

      return instance_params;
    }
  }
  
};

template<typename Graph,
	 typename work_item_t,
	 typename MessageGenerator = 
	 amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
struct triangle_verifier {

  typedef triangle_verifier<Graph, MessageGenerator> self_type;
  typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  struct verify_handler;

  struct minimum_pair_first
  {
    template<typename T>
    const T& operator()(const T& x, const T& y) const { return x < y ? x : y; }

    template<typename F>
    struct result {
      typedef typename boost::function_traits<F>::arg1_type type;
    };
  };


  typedef typename MessageGenerator::template call_result<work_item_t, 
							  verify_handler, 
							  owner_from_pair<OwnerMap, work_item_t>, 
							  amplusplus::idempotent_combination_t<minimum_pair_first > >::type RelaxMessage;


public:
  triangle_verifier(Graph& gr,
		    std::vector<work_item_t>& t,
		    amplusplus::transport& tr,
		    int offset,
		    MessageGenerator message_gen =
		    MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12))):
                                dummy_first_member_for_init_order
				((amplusplus::register_mpi_datatype<work_item_t>(), 
				  0)),
				g(gr),
				owner(boost::get(boost::vertex_owner, g)), 
				triangles(t),
				transport(tr),
				core_offset(offset),
				relax_msg(message_gen, 
					  transport, 
					  owner_from_pair<OwnerMap, work_item_t>(owner), 
					  amplusplus::idempotent_combination(minimum_pair_first())){}
    

  void initialize() {
    relax_msg.set_handler(verify_handler(*this, g));
  }

  void operator()(int tid) {
    AMPLUSPLUS_WITH_THREAD_ID(tid) {
      int nthreads = transport.get_nthreads();
      if (0 == tid) {
	// Set the number of threads to the barrier
	t_bar.reset(new amplusplus::detail::barrier(nthreads));
      }

      { amplusplus::scoped_epoch epoch(transport); }

      // Now above if branch needs to be executed to every thread
      // Therefore wait till every thread comes to this point
      t_bar->wait();

      // if two processes are running on the same node, core_offset
      // is important to achieve thread affinity
      if (pin(tid+core_offset) != 0) {
	std::cerr << "[ERROR] Unable to pin current thread to "
		  << "core : " << tid << std::endl;
	assert(false);
      }

      // wait till all threads are pinned
      t_bar->wait();
      { amplusplus::scoped_epoch epoch(transport); }

      validate_thread_core_relation();

      {
	amplusplus::scoped_epoch epoch(transport);

	for(int i=tid; i < triangles.size(); i=(i+nthreads)) {
	  work_item_t wi = triangles[i];

	  Vertex succ = wi.first;
	  Vertex cur = wi.second.first;
	  Vertex pred = wi.second.second;

	  assert(boost::get(owner, succ) == transport.rank());

	  bool fedgecur = false;
	  bool fedgepred = false;

	  //check whether there are edges from succ to cur and pred
	  BGL_FORALL_OUTEDGES_T(succ, e, g, Graph) {
	    Vertex u = boost::target(e, g);
	    if (u == cur) {
	      fedgecur = true;
	      if (fedgepred)
		break;
	    }

	    if (u == pred) {
	      fedgepred = true;
	      if (fedgecur)
		break;
	    }
	  }

	  if (fedgecur && fedgepred) {
	    // we know there is an edge from succ to cur and pred
	    // we now need to find whether there is an edge between pred and cur
	    // if one of pred and cur in current locality we dont need distributed
	    // communication. If not we need to send a message to owner of cur or pred
	    bool foundlocal = false;
	    bool remote = false;
	    if (boost::get(owner, cur) == transport.rank()) {
	      BGL_FORALL_OUTEDGES_T(cur, e, g, Graph) {
		Vertex u = boost::target(e, g);
		if (u == pred) {
		  foundlocal = true;
		  break;
		}
	      }
	    } else if (boost::get(owner, pred) == transport.rank()) {
	      BGL_FORALL_OUTEDGES_T(pred, e, g, Graph) {
		Vertex u = boost::target(e, g);
		if (u == cur) {
		  foundlocal = true;
		  break;
		}
	      }
	    } else {
	      // send message
	      remote = true;

	      work_item_t wisend;
	      wisend.first = wi.second.first;
	      wisend.second.first = wi.second.second;
	      wisend.second.second = wi.first;
	  
	      relax_msg.send(wisend);
	    }
	
	    if ((!remote) && (!foundlocal)) {
	      std::cout << "[ERROR] Did not find an edge between vertices, ("
			<< cur << ", " << pred << "), "
			<< "for the triangle (" << succ << ", "
			<< cur << ", " << pred << ")" << std::endl;
	      assert(false);
	    }
	  } else {
	    if (!fedgecur) {
	      std::cout << "[ERROR] Did not find an edge between vertices, ("
			<< succ << ", " << cur << "), "
			<< "for the triangle (" << succ << ", "
			<< cur << ", " << pred << ")" << std::endl;
	    }

	    if (!fedgepred) {
	      std::cout << "[ERROR] Did not find an edge between vertices, ("
			<< succ << ", " << pred << "), "
			<< "for the triangle (" << succ << ", "
			<< cur << ", " << pred << ")" << std::endl;
	    }

	    assert(false);
	  }
	} // end for loop
      } // end epoch
    }// end of ampp thread id
  } // end of operator

private:
  int dummy_first_member_for_init_order;
  std::vector<work_item_t>& triangles;
  amplusplus::transport& transport;
  const Graph& g;
  const OwnerMap& owner;
  int core_offset;
  RelaxMessage relax_msg;
  boost::shared_ptr<amplusplus::detail::barrier> t_bar;

public:
  struct verify_handler {
    verify_handler() : self(NULL), g(NULL){}
    verify_handler(triangle_verifier& slf,
		   Graph& gi) : self(&slf),
			       g(&gi){}

    void operator()(const work_item_t& data) const {
      Vertex cur = data.first;
      Vertex pred = data.second.first;
      Vertex succ = data.second.second;

      bool foundedge = false;
      // is there an edge between cur and pred, then we are good
      BGL_FORALL_OUTEDGES_T(cur, e, *g, Graph) {
	Vertex u = boost::target(e, *g);
	if (u == pred) {
	  foundedge = true;
	  break;
	}
      }

      if (!foundedge) {
	std::cout << "[ERROR] Did not find an edge between vertices, ("
		  << cur << ", " << pred << "), "
		  << "for the triangle (" << succ << ", "
		  << cur << ", " << pred << ")" << std::endl;
	assert(false);
      }
    }

  protected:
    triangle_verifier* self;
    const Graph* g;
  };
};


class TCExecutor {
private:

  template <typename Graph, typename work_item_t>
  bool verify_tc(amplusplus::transport& trans,  
		 const instance_params& runtime_params,
		 Graph& g,
		 std::vector<work_item_t>& triangles) {

    time_type start_time = get_time();

    trans.set_nthreads(runtime_params.threads);
    triangle_verifier<Graph, work_item_t> tv(g, triangles, trans, sched_getcpu());
    tv.initialize();

    //for each work item check whether there is an edge between them
    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
    for (int i = 0; i < runtime_params.threads - 1; ++i) {
      boost::thread thr(boost::ref(tv), i + 1);
      threads[i].swap(thr);
    }
	  
    tv(0);

    for (int i = 0; i < runtime_params.threads - 1; ++i)
      threads[i].join();

    trans.set_nthreads(1);

    time_type end_time = get_time();

    std::cout << "[INFO] Verification took : " << (end_time - start_time)
	      << std::endl;

    return true;
  }

  template<typename graph_create_params>
  void simple_count_verify(uint64_t triangles,
			   graph_create_params& gparams,
			   tc_instance_params& tc_params) {
    if (gparams.preprocess) {
      if (_RANK == 0) {
	std::cout << "Expected triangles : " << tc_params.triangles
		  << ", counted : " << triangles
		  << std::endl;
	assert(triangles == tc_params.triangles);
      }
    }
  }

  // Triangle Counting algorithms
  template <typename Graph, 
	    typename IdDistribution, 
	    typename MessageGenerator,
	    typename graph_create_params>
  time_type
  run_tc(amplusplus::transport& trans, 
	 amplusplus::transport& barrier_trans, 
	 Graph& g,  
	 const IdDistribution& idd,
	 MessageGenerator msg_gen, 
	 graph_create_params& gparams,
	 instance_params& runtime_params,
	 tc_instance_params& tc_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;


    if (trans.rank() == 0)
      std::cout << "Initializing triangle counting map ..." << std::endl;

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::triangle_counting<Graph, 
						 IdDistribution,
						 MessageGenerator>
      D(g, trans, idd, sched_getcpu(), msg_gen);
    
    trans.set_nthreads(1);

    { amplusplus::scoped_epoch epoch(barrier_trans); }

    // Many threads now
    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Invoking algorithm ..." << std::endl;

    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
    for (int i = 0; i < runtime_params.threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < runtime_params.threads - 1; ++i)
      threads[i].join();
	  
    time_type end = get_time();

    if (trans.rank() == 0)
      std::cout << "Algorithm done ..." << std::endl;

    time_type start = D.get_start_time();

    D.print_triangle_counts();

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

#ifdef TRIANGLE_ENUMERATE
    typedef typename boost::graph::distributed::triangle_counting<Graph, 
						     IdDistribution,
						     MessageGenerator>::work_item_t work_item_t;

    std::vector<work_item_t> local_triangles;
    D.get_local_triangles(local_triangles);
#endif

    // Calculate number of triangles
    uint64_t globaltraingles = 0;
    uint64_t localtriangles = D.get_local_triangle_counts();
    // MPI reduce to find the all triangles
    MPI_Reduce(&localtriangles, &globaltraingles, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);


    if (tc_params.verify) {
      if (trans.rank()==0)
	std::cout << "Verifying counts ..." << std::endl;

      simple_count_verify(globaltraingles, gparams, tc_params);

#ifdef TRIANGLE_ENUMERATE
      if (!verify_tc<Graph, work_item_t>(trans, runtime_params, 
					 g, local_triangles)) {
	std::cout << "Triangle counting Verification Failed" << std::endl;
	assert(false);
	return 0;
      }
#endif

    }

    if (trans.rank() == 0)
      std::cout << "Total number of triangles counted : " << globaltraingles
		<< ", in " << (end-start) << " time." << std::endl;
	  
    return D.get_elapsed_time();
  }


  // Triangle Counting algorithms
  template <typename Graph, 
	    typename IdDistribution, 
	    typename MessageGenerator,
	    typename graph_create_params>
  time_type
  run_tc_blocked(amplusplus::transport& trans, 
		 amplusplus::transport& barrier_trans, 
		 Graph& g,  
		 const IdDistribution& idd,
		 MessageGenerator msg_gen, 
		 graph_create_params& gparams,
		 instance_params& runtime_params,
		 tc_instance_params& tc_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

    typedef std::vector<Vertex> InternalContainerType;
    typedef std::vector<InternalContainerType> NeighborsType;
    typedef boost::iterator_property_map<typename NeighborsType::iterator, VertexIndexMap>  NeighborsMap;

    NeighborsType vec_predecessors(num_vertices(g), InternalContainerType(0));
    NeighborsMap predecessor_map(vec_predecessors.begin(), get(boost::vertex_index, g));

    NeighborsType vec_successors(num_vertices(g), InternalContainerType(0));
    NeighborsMap successor_map(vec_successors.begin(), get(boost::vertex_index, g));

    if (trans.rank() == 0)
      std::cout << "[INFO] Initializing triangle counting map ..." << std::endl;

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "[INFO] Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::triangle_counting_blocked<Graph, 
							 IdDistribution,
							 NeighborsMap>
      D(g, trans, idd, 
	sched_getcpu(), 
	tc_params.block_size, 
	tc_params.suc_block_size, 
	runtime_params.coalescing, 
	predecessor_map,
	successor_map);
    
    trans.set_nthreads(1);

    { amplusplus::scoped_epoch epoch(barrier_trans); }

    // Many threads now
    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "[INFO] Invoking algorithm ..." << std::endl;

    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
    for (int i = 0; i < runtime_params.threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < runtime_params.threads - 1; ++i)
      threads[i].join();
	  
    time_type end = get_time();

    if (trans.rank() == 0)
      std::cout << "[INFO] Algorithm done ..." << std::endl;

    time_type start = D.get_start_time();

    //D.print_triangle_counts();
    //std::cout << trans.rank() << " Done printing counts" << std::endl;

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

#ifdef TRIANGLE_ENUMERATE
    typedef typename boost::graph::distributed::triangle_counting<Graph, 
						     IdDistribution,
						     MessageGenerator>::work_item_t work_item_t;

    std::vector<work_item_t> local_triangles;
    D.get_local_triangles(local_triangles);
#endif

    // Calculate number of triangles
    uint64_t globaltraingles = 0;
    uint64_t localtriangles = D.get_local_triangle_counts();
    // MPI reduce to find the all triangles
    MPI_Reduce(&localtriangles, &globaltraingles, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);

    if (trans.rank()==0) {
      if (tc_params.verify) {
	std::cout << "Verifying counts ..." << std::endl;

	simple_count_verify(globaltraingles, gparams, tc_params);

#ifdef TRIANGLE_ENUMERATE
	if (!verify_tc<Graph, work_item_t>(trans, runtime_params, 
					   g, local_triangles)) {
	  std::cout << "Triangle counting Verification Failed" << std::endl;
	  assert(false);
	  return 0;
	}
#endif

      }

      std::cout << "Total number of triangles counted : " << globaltraingles
		<< ", in " << (end-start) << " time." << std::endl;
    }

#ifdef TC_STATS
    D.print_stats();
#endif
	  
    return D.get_elapsed_time();
  }

  // Triangle Counting algorithms
  template <typename Graph, 
	    typename IdDistribution, 
	    typename MessageGenerator,
	    typename graph_create_params>
  time_type
  run_tc_sps(amplusplus::transport& trans, 
	     amplusplus::transport& barrier_trans, 
	     Graph& g,  
	     const IdDistribution& idd,
	     MessageGenerator msg_gen, 
	     graph_create_params& gparams,
	     instance_params& runtime_params,
	     tc_instance_params& tc_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

    typedef std::vector<Vertex> InternalContainerType;
    typedef std::vector<InternalContainerType> NeighborsType;
    typedef boost::iterator_property_map<typename NeighborsType::iterator, VertexIndexMap>  NeighborsMap;

    NeighborsType vec_predecessors(num_vertices(g), InternalContainerType(0));
    NeighborsMap predecessor_map(vec_predecessors.begin(), get(boost::vertex_index, g));

    NeighborsType vec_successors(num_vertices(g), InternalContainerType(0));
    NeighborsMap successor_map(vec_successors.begin(), get(boost::vertex_index, g));

    if (trans.rank() == 0)
      std::cout << "[INFO] Initializing triangle counting map ..." << std::endl;

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "[INFO] Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::triangle_counting_sps<Graph, 
							 IdDistribution,
							 NeighborsMap>
      D(g, trans, idd, 
	sched_getcpu(), 
	tc_params.block_size, 
	tc_params.suc_block_size, 
	runtime_params.coalescing, 
	predecessor_map,
	successor_map);
    
    trans.set_nthreads(1);

    { amplusplus::scoped_epoch epoch(barrier_trans); }

    // Many threads now
    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "[INFO] Invoking algorithm ..." << std::endl;

    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
    for (int i = 0; i < runtime_params.threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < runtime_params.threads - 1; ++i)
      threads[i].join();
	  
    time_type end = get_time();

    if (trans.rank() == 0)
      std::cout << "[INFO] Algorithm done ..." << std::endl;

    time_type start = D.get_start_time();

    //D.print_triangle_counts();
    //std::cout << trans.rank() << " Done printing counts" << std::endl;

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

#ifdef TRIANGLE_ENUMERATE
    typedef typename boost::graph::distributed::triangle_counting<Graph, 
						     IdDistribution,
						     MessageGenerator>::work_item_t work_item_t;

    std::vector<work_item_t> local_triangles;
    D.get_local_triangles(local_triangles);
#endif

    // Calculate number of triangles
    uint64_t globaltraingles = 0;
    uint64_t localtriangles = D.get_local_triangle_counts();
    // MPI reduce to find the all triangles
    MPI_Reduce(&localtriangles, &globaltraingles, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);

    if (trans.rank()==0) {
      if (tc_params.verify) {
	std::cout << "Verifying counts ..." << std::endl;

	simple_count_verify(globaltraingles, gparams, tc_params);

#ifdef TRIANGLE_ENUMERATE
	if (!verify_tc<Graph, work_item_t>(trans, runtime_params, 
					   g, local_triangles)) {
	  std::cout << "Triangle counting Verification Failed" << std::endl;
	  assert(false);
	  return 0;
	}
#endif

      }

      std::cout << "Total number of triangles counted : " << globaltraingles
		<< ", in " << (end-start) << " time." << std::endl;
    }

#ifdef TC_STATS
    D.print_stats();
#endif
	  
    return D.get_elapsed_time();
  }


  // Triangle Counting algorithms
  template <typename Graph, 
	    typename IdDistribution, 
	    typename MessageGenerator,
	    typename graph_create_params>
  time_type
  run_tc_sucsuc(amplusplus::transport& trans, 
		 amplusplus::transport& barrier_trans, 
		 Graph& g,  
		 const IdDistribution& idd,
		 MessageGenerator msg_gen, 
		 graph_create_params& gparams,
		 instance_params& runtime_params,
		 tc_instance_params& tc_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;

    typedef std::vector<Vertex> InternalContainerType;
    typedef std::vector<InternalContainerType> NeighborsType;
    typedef boost::iterator_property_map<typename NeighborsType::iterator, VertexIndexMap>  NeighborsMap;

    NeighborsType vec_successors(num_vertices(g), InternalContainerType(0));
    NeighborsMap successor_map(vec_successors.begin(), get(boost::vertex_index, g));

    if (trans.rank() == 0)
      std::cout << "Initializing triangle counting map ..." << std::endl;

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::triangle_counting_sucsuc<Graph, 
							 IdDistribution,
							 NeighborsMap>
      D(g, trans, idd, 
	sched_getcpu(), 
	tc_params.block_size, 
	tc_params.suc_block_size, 
	runtime_params.coalescing, 
	successor_map);
    
    trans.set_nthreads(1);

    { amplusplus::scoped_epoch epoch(barrier_trans); }

    // Many threads now
    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Invoking algorithm ..." << std::endl;

    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
    for (int i = 0; i < runtime_params.threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < runtime_params.threads - 1; ++i)
      threads[i].join();
	  
    time_type end = get_time();

    if (trans.rank() == 0)
      std::cout << "Algorithm done ..." << std::endl;

    time_type start = D.get_start_time();

    D.print_triangle_counts();
    std::cout << trans.rank() << " Done printing counts" << std::endl;

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

#ifdef TRIANGLE_ENUMERATE
    typedef typename boost::graph::distributed::triangle_counting<Graph, 
						     IdDistribution,
						     MessageGenerator>::work_item_t work_item_t;

    std::vector<work_item_t> local_triangles;
    D.get_local_triangles(local_triangles);
#endif

    // Calculate number of triangles
    uint64_t globaltraingles = 0;
    uint64_t localtriangles = D.get_local_triangle_counts();
    // MPI reduce to find the all triangles
    MPI_Reduce(&localtriangles, &globaltraingles, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);

    if (trans.rank()==0) {
      if (tc_params.verify) {
	std::cout << "Verifying counts ..." << std::endl;

	simple_count_verify(globaltraingles, gparams, tc_params);

#ifdef TRIANGLE_ENUMERATE
	if (!verify_tc<Graph, work_item_t>(trans, runtime_params, 
					   g, local_triangles)) {
	  std::cout << "Triangle counting Verification Failed" << std::endl;
	  assert(false);
	  return 0;
	}
#endif

      }

      std::cout << "Total number of triangles counted : " << globaltraingles
		<< ", in " << (end-start) << " time." << std::endl;
    }

#ifdef TC_STATS
    D.print_stats();
#endif
	  
    return D.get_elapsed_time();
  }


  // Triangle Counting algorithms
  template <typename Graph, 
	    typename IdDistribution, 
	    typename MessageGenerator,
	    typename graph_create_params>
  time_type
  run_tc_level(amplusplus::transport& trans, 
		 amplusplus::transport& barrier_trans, 
		 Graph& g,  
		 const IdDistribution& idd,
		 MessageGenerator msg_gen, 
		 graph_create_params& gparams,
		 instance_params& runtime_params,
		 tc_instance_params& tc_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    typedef boost::iterator_property_map<typename std::vector<int>::iterator, VertexIndexMap> FlagMap;
    std::vector<int> boolvec(num_vertices(g), 0);
    FlagMap level1flags(boolvec.begin(), get(boost::vertex_index, g));

    if (trans.rank() == 0)
      std::cout << "Initializing triangle counting map ..." << std::endl;

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::triangle_counting_level<Graph, FlagMap,
							 IdDistribution,
							 MessageGenerator>
      D(g, trans, idd, level1flags, sched_getcpu(), msg_gen);
    
    trans.set_nthreads(1);

    { amplusplus::scoped_epoch epoch(barrier_trans); }

    // Many threads now
    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Invoking algorithm ..." << std::endl;

    boost::scoped_array<boost::thread> threads(new boost::thread[runtime_params.threads - 1]);
    for (int i = 0; i < runtime_params.threads - 1; ++i) {
      boost::thread thr(boost::ref(D), i + 1);
      threads[i].swap(thr);
    }
	  
    D(0);
    
    for (int i = 0; i < runtime_params.threads - 1; ++i)
      threads[i].join();
	  
    time_type end = get_time();

    if (trans.rank() == 0)
      std::cout << "Algorithm done ..." << std::endl;

    time_type start = D.get_start_time();

    D.print_triangle_counts();
    std::cout << trans.rank() << " Done printing counts" << std::endl;

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

    // Calculate number of triangles
    uint64_t globaltraingles = 0;
    uint64_t localtriangles = D.get_local_triangle_counts();
    // MPI reduce to find the all triangles
    MPI_Reduce(&localtriangles, &globaltraingles, 
	       1, MPI_LONG_LONG_INT , MPI_SUM, 0, MPI_COMM_WORLD);

    if (trans.rank()==0) {
      if (tc_params.verify) {
	std::cout << "Verifying counts ..." << std::endl;
	simple_count_verify(globaltraingles, gparams, tc_params);
      }

      std::cout << "Total number of triangles counted : " << globaltraingles
		<< ", in " << (end-start) << " time." << std::endl;
    }
	  
    return D.get_elapsed_time();
  }



public:
  template <typename Graph, typename MessageGenerator, typename graph_create_params>
  time_type operator()(const Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       tc_instance_params& tc_params) { // remove

    if (tc_params.algorithm == dag_tc_striped) {
      amplusplus::transport barrier_trans = trans.clone();
      if (tc_params.id_distribution == vertical) {
	return run_tc(trans,
		      barrier_trans,
		      g,
		      block_id_distribution<Graph>(g, gparams.n),
		      msg_gen,
		      gparams,
		      runtime_params,
		      tc_params); 
      } else if (tc_params.id_distribution == horizontal) {
	return run_tc(trans,
		      barrier_trans,
		      g,
		      row_id_distribution<Graph>(g, trans.size()),
		      msg_gen,
		      gparams,
		      runtime_params,
		      tc_params); 

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
    } else if (tc_params.algorithm == dag_tc_blocked) {
      amplusplus::transport barrier_trans = trans.clone();
      if (tc_params.id_distribution == vertical) {
	return run_tc_blocked(trans,
			      barrier_trans,
			      g,
			      block_id_distribution<Graph>(g, gparams.n),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 
      } else if (tc_params.id_distribution == horizontal) {
	return run_tc_blocked(trans,
			      barrier_trans,
			      g,
			      row_id_distribution<Graph>(g, trans.size()),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 

      } else if (tc_params.id_distribution == degree) {
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	return run_tc_blocked(trans,
			      barrier_trans,
			      g,
			      degree_id_distribution<Graph>(g, trans.size(), tc_params.degreeasc),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
    } else if (tc_params.algorithm == dag_tc_sps) {
      amplusplus::transport barrier_trans = trans.clone();
      if (tc_params.id_distribution == vertical) {
	return run_tc_sps(trans,
			  barrier_trans,
			  g,
			  block_id_distribution<Graph>(g, gparams.n),
			  msg_gen,
			  gparams,
			  runtime_params,
			  tc_params); 
      } else if (tc_params.id_distribution == horizontal) {
	return run_tc_sps(trans,
			  barrier_trans,
			  g,
			  row_id_distribution<Graph>(g, trans.size()),
			  msg_gen,
			  gparams,
			  runtime_params,
			  tc_params); 

      } else if (tc_params.id_distribution == degree) {
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	return run_tc_sps(trans,
			  barrier_trans,
			  g,
			  degree_id_distribution<Graph>(g, trans.size(), tc_params.degreeasc),
			  msg_gen,
			  gparams,
			  runtime_params,
			  tc_params); 

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
    } else if (tc_params.algorithm == dag_tc_sucsuc) {
      amplusplus::transport barrier_trans = trans.clone();
      if (tc_params.id_distribution == vertical) {
	return run_tc_sucsuc(trans,
			      barrier_trans,
			      g,
			      block_id_distribution<Graph>(g, gparams.n),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 
      } else if (tc_params.id_distribution == horizontal) {
	return run_tc_sucsuc(trans,
			      barrier_trans,
			      g,
			      row_id_distribution<Graph>(g, trans.size()),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 

      } else if (tc_params.id_distribution == degree) {
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	return run_tc_sucsuc(trans,
			      barrier_trans,
			      g,
			      degree_id_distribution<Graph>(g, trans.size(), tc_params.degreeasc),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
    } else if (tc_params.algorithm == dag_tc_level) {
      amplusplus::transport barrier_trans = trans.clone();
      if (tc_params.id_distribution == vertical) {
	return run_tc_level(trans,
			      barrier_trans,
			      g,
			      block_id_distribution<Graph>(g, gparams.n),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 
      } else if (tc_params.id_distribution == horizontal) {
	return run_tc_level(trans,
			      barrier_trans,
			      g,
			      row_id_distribution<Graph>(g, trans.size()),
			      msg_gen,
			      gparams,
			      runtime_params,
			      tc_params); 

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
    } else {
      std::cerr << "Only supports triangle counting synchronous version (--run_tc_striped or --run_tc_blocked)"
		<< std::endl;
      assert(false);
    }
  }
};


int main(int argc, char* argv[]) {
  //std::cout << "printing core id for process ..." << std::endl;
  print_core_id();

  /*{
    std::cout << "Allocating memory in main" << std::endl;
    typedef std::pair<uint64_t, uint64_t> edge_t; 
    typedef std::vector<edge_t> edges_t;
    edges_t edges;
    int64_t directedes = 1468365182;
    int64_t total_alloc = directedes * (int64_t)2; 
    edges.reserve(total_alloc);
    std::cout << "Allocated memory in main" << std::endl;
    }*/

  executor<TCExecutor, tc_params, tc_instance_params> tc_executor;
  tc_executor.execute(argc, argv);  
}
