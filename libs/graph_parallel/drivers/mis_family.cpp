// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine

//======== Maximal Independent Set Algortihms================//
// Driver for MIS family of algorithms.
//===========================================================//


#include <iostream>
#include "common/synthetic_generator.hpp"
#include "common/parser.hpp"
#include "common/executor.hpp"
#include "common/vertex_permutations.hpp"
#include <boost/graph/distributed/mis_delta.hpp>
#include <boost/graph/distributed/mis.hpp>
#include <boost/graph/distributed/luby_mis.hpp>
#include <boost/graph/distributed/thread_pq_def.hpp>


enum mis_algorithm {
  fix,
  fix_bucket,
  luby_a,
  luby_av1,
  luby_av2,
  luby_b
};


class mis_instance_params {

public:
  mis_algorithm algorithm;
  id_distribution_t id_distribution;
  bool verify;
  mis_instance_params(mis_algorithm& alg,
		      id_distribution_t& idd,
		      bool v):
    algorithm(alg), id_distribution(idd), verify(v){}

  void print() {
      std::cout << "Algorithm -- ";
      if (algorithm == fix) 
	std::cout << "FIX";
      else if(algorithm == fix_bucket)
	std::cout << "FIX-Bucket";
      else if(algorithm == luby_a)
	std::cout << "Luby A";
      else if(algorithm == luby_av1)
	std::cout << "Luby AV1";
      else if(algorithm == luby_av2)
	std::cout << "Luby AV2";
      else if(algorithm == luby_b)
	std::cout << "Luby B";
      else
	std::cerr << "Invalid !" << std::endl;

      std::cout << std::endl;


    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "verify : " << verify << std::endl;
  }

  std::string get_algorithm() {
    if (algorithm == fix) 
      return "FIX";
    else if(algorithm == fix_bucket)
      return "FIX-Bucket";
    else if(algorithm == luby_a)
      return "Luby A";
    else if(algorithm == luby_av1)
      return "Luby AV1";
    else if(algorithm == luby_av2)
      return "Luby AV2";
    else if(algorithm == luby_b)
      return "Luby B";
    else
      return "Invalid !";
  }

};


class mis_params {

private:
  std::vector<mis_instance_params> instance_params;
  std::vector<mis_algorithm> algorithms;
  id_distribution_t id_distribution = horizontal; // default
  bool verify = false;

public:
  bool parse(int argc, char* argv[]){
    for (int i = 1; i < argc; ++i) {
      mis_algorithm algorithm;
      std::string arg = argv[i];
      if (arg == "--run_fix_mis") {
	algorithm = fix;
	algorithms.push_back(algorithm);
      }

      if (arg == "--run_bucket_mis") {
	algorithm = fix_bucket;
	algorithms.push_back(algorithm);
      }

      if (arg == "--id-distribution") {
	if (strcmp(argv[i+1],"vertical") == 0)
	  id_distribution = vertical;
	else if (strcmp(argv[i+1],"horizontal") == 0)
	  id_distribution = horizontal;
	else {
	  std::cout << "Invalid id distribution type. Available types are vertical and horizontal" 
		    << std::endl;
	  return false;
	}
      }

      if (arg == "--verify") {
	verify = true;
      }

      if (arg == "--luby_algorithms") {
	std::vector<std::string> luby_algorithms;
	luby_algorithms = extract_params<std::string> ( argv[i+1] );
	BOOST_FOREACH(std::string al, luby_algorithms) {
	  if (al == "A") 
	    algorithm = luby_a;
	  else if (al == "AV1")
	    algorithm = luby_av1;
	  else if (al == "AV2")
	    algorithm = luby_av2;
	  else if (al == "B")
	    algorithm = luby_b;
	  else {
	    std::cerr << "Invalid Luby algorithm. Available algorithms --"
		      << "A, AV1, AV2, B" << std::endl;
	    return false;
	  }

	  algorithms.push_back(algorithm);
	}
      }
    }
  }

  void print() {
    BOOST_FOREACH(mis_algorithm alg, algorithms) {
      std::cout << "Algorithm -- ";
      if (alg == fix) 
	std::cout << "FIX";
      else if(alg == fix_bucket)
	std::cout << "FIX-Bucket";
      else if(alg == luby_a)
	std::cout << "Luby A";
      else if(alg == luby_av1)
	std::cout << "Luby AV1";
      else if(alg == luby_av2)
	std::cout << "Luby AV2";
      else if(alg == luby_b)
	std::cout << "Luby B";
      else
	std::cerr << "Invalid !" << std::endl;

      std::cout << std::endl;
    }

    if (id_distribution == vertical)
      std::cout << "id distribution : vertical" << std::endl;

    if (id_distribution == horizontal)
      std::cout << "id distribution : horizontal" << std::endl;
    
    std::cout << "verify : " << verify << std::endl;
  }

  const std::vector<mis_instance_params>&
  get_instance_params() {
    if (instance_params.empty()) {
      BOOST_FOREACH(mis_algorithm alg, algorithms) {
	instance_params.push_back(mis_instance_params(alg, id_distribution, verify));
      }
    }

    return instance_params;
  }
  
};

class MISExecutor {
private:

  template <typename Graph, typename MISMap>
  bool verify_mis(amplusplus::transport& trans,  Graph& g, MISMap& mis) {
    typedef typename boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
    OwnerMap owner(get(boost::vertex_owner, g));

    mis.set_consistency_model(boost::parallel::cm_forward || boost::parallel::cm_backward);
    mis.set_max_ghost_cells(0);

    if (trans.rank() == 0) std::cout<<"Verifying MIS results......";

    {
      amplusplus::scoped_epoch epoch(g.transport());
	  
      BGL_FORALL_VERTICES_T(v, g, Graph) {
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  get(mis, target(e, g));
	}
      }
    }
	    
    bool result = true;
    int found_mis = 0;
#ifdef PRINT_DEBUG
    std::cout << "MIS = {";
#endif
    BGL_FORALL_VERTICES_T(v, g, Graph) {
      if (get(mis, v) ==  MIS_FIX1) {// v in mis, none of the neigbours should be in mis
#ifdef PRINT_DEBUG
	if (v == 35) {
	  std::cout << "Printing all neighbors of " << v
		    << " -- [";

	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    std::cout << target(e, g) << ", ";
	  }
	  std::cout << "]" << std::endl;
	}

	std::cout << v << "-[";
#endif
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
#ifdef PRINT_DEBUG
	  std::cout << target(e, g) << ", ";
#endif
	  if (v == target(e, g))
	    continue;

	  if (get(mis, target(e, g)) == MIS_FIX1) {
	    std::cout << "[FAIL] Vertex : " << v << " and its neigbour " << target(e, g)
		      << " are both in MIS." 
		      << " Vertex value : " << get(mis, v)
		      << " neighbor value : " << get(mis, target(e, g)) 
		      << std::endl;
	    std::cout << "Neighbors of " << v << "- {";
	    BGL_FORALL_ADJ_T(v, u1, g, Graph) {
	      std::cout << "(" << u1 << ", " << get(mis, u1) << ")";
	    }
	    std::cout << "}" << std::endl;

	    auto k = target(e, g);
	    if (get(owner, k) == g.transport().rank()) {
	      std::cout << "Neighbors of " << k << "- {";
	      BGL_FORALL_ADJ_T(k, u2, g, Graph) {
		std::cout << "(" << u2 << ", " << get(mis, u2) << ")";
	      }
	      std::cout << "}" << std::endl;
	    }

	    result = false;
	  } else {
#ifdef PRINT_DEBUG
	    if (get(mis, target(e, g)) != MIS_FIX0) {
	      std::cout << "vertex : " 
			<< target(e, g) 
			<< ", is in " 
			<< get(mis, target(e, g))
			<< std::endl;
	    }
#endif

	    // cannot be MIS_UNFIX
	    assert(get(mis, target(e, g)) == MIS_FIX0);
	    found_mis = 1;
	  }
	}
#ifdef PRINT_DEBUG
	std::cout << "], ";
#endif
      } else {
	// if v is not in mis, at least one of its neighbors must
	// be in mis
	if (get(mis, v) != MIS_FIX0) {
	  std::cout << "[ERROR] Vertex - " << v << " is in " << get(mis, v)
		    << std::endl;
	}

	if (get(mis, v) != MIS_FIX0) {
	  std::cout << "Error : " << v << " is not in MIS_FIX0. But in " << get(mis, v)
		    << ", owner of " << v << ":" << get(owner, v)
		    << ", current rank : " << _RANK << std::endl;
	  BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	    auto k = target(e, g);
	    std::cout << "neigbor : " << k << " mis : " << get(mis, k)
		      << std::endl;
	  }
	}
	assert(get(mis, v) == MIS_FIX0);

	bool inmis = false;
	BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	  if (v == target(e, g))
	    continue;

	  if (get(mis, target(e, g)) == MIS_FIX1) {
	    inmis = true;
	    break;
	  }
	}
      
	if (!inmis) {
	  std::cout << "Vertex : " << v << " and none of its neighbors is in MIS" << std::endl;
	  assert(false);
	}
      }
    }
#ifdef PRINT_DEBUG
    std::cout << "}" << std::endl;
#endif

    // Run a reduction on found mis.
    // We cannot expect every rank will have found_mis true.
    // E.g:- rank0 - {0,1}, rank1 - {281474976710656, 281474976710657}
    // If every vertex is connect to each other, a vertex (Lets say 0) will
    // be in one rank. In the second rank all vertices are out of mis.
    int or_found_mis = 0;
    MPI_Allreduce(&found_mis, &or_found_mis, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (or_found_mis == 0) {
      std::cout << "Did not find any MIS" <<std::endl;
      result = false;
    }

    return result;
  }


  // Luby Algorithms
  template <typename SelectGenerator, typename Graph, 
	    typename graph_create_params,
	    typename MessageGenerator>
  time_type
  run_luby_maximal_is(amplusplus::transport& trans, 
		    amplusplus::transport& barrier_trans, 
		    Graph& g,  
		    MessageGenerator msg_gen, 
		    graph_create_params& gparams,
		    instance_params& runtime_params,
		    mis_instance_params& mis_params) {

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;

    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    typedef boost::iterator_property_map<typename std::vector<random_t>::iterator, VertexIndexMap>  RandomMap;

    std::vector<state_t> misvec(num_vertices(g), MIS_UNFIX);
    typedef boost::iterator_property_map<typename std::vector<state_t>::iterator, VertexIndexMap>  MISMap;
    MISMap mis(misvec.begin(), get(boost::vertex_index, g));


    // create a property map
    std::vector<random_t> pivec(num_vertices(g), 0);
    RandomMap rmap(pivec.begin(), get(boost::vertex_index, g)); // TODO remove copying

    if (trans.rank() == 0)
      std::cout << "Initializing mis map ..." << std::endl;

    BGL_FORALL_VERTICES_T(v, g, Graph) 
      { put(mis, v, MIS_UNFIX); }

    // TODO find why we need this ?
    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::luby_mis<Graph, MISMap, RandomMap, 
					SelectGenerator,
					MessageGenerator>
      D(g, mis, rmap, gparams.n, trans, sched_getcpu(), msg_gen);
    
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

#ifdef MIS_STATS
    D.print_stats();
#endif

    time_type start = D.get_start_time();

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

    if (mis_params.verify) {
      if (trans.rank()==0)
	std::cout << "Verifying mis ..." << std::endl;

      if (!verify_mis(trans, g, mis)) {
	std::cout << "MIS Verification Failed" << std::endl;
	assert(false);
	return 0;
      }
    }

    vertices_size_type visited = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (get(mis, v) != MIS_UNFIX)
	++visited; 
    }
	  
    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
      r(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r(visited);
	  
    if (mis_params.verify)
      if (trans.rank() == 0)
	std::cout << "Visited " << total << " vertices of " << gparams.n << " in " << print_time(end - start) 
		  << std::endl;

	  
    return end - start;
  }


  // FIX algorithms
  template <typename Graph, typename IdDistribution, typename MessageGenerator, 
	    typename graph_create_params,
	    typename PriorityQueueGenerator = boost::graph::distributed::thread_priority_queue_gen>
  time_type
  run_fix_mis(amplusplus::transport& trans, 
	      amplusplus::transport& barrier_trans, 
	      Graph& g,  
	      const IdDistribution& idd,
	      MessageGenerator msg_gen, 
	      graph_create_params& gparams,
	      instance_params& runtime_params,
	      mis_instance_params& mis_params) { 

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;

#ifdef MIS_DISTANCE_ORDERING
    std::cout << "Running with distance ordering ..." << std::endl;
#endif

    if (trans.rank() == 0)
      std::cout << "Initializing mis map ..." << std::endl;

    std::vector<state_t> misvec(num_vertices(g), MIS_UNFIX);
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    typedef boost::iterator_property_map<typename std::vector<state_t>::iterator, VertexIndexMap>  MISMap;
    MISMap mis(misvec.begin(), get(boost::vertex_index, g));

    BGL_FORALL_VERTICES_T(v, g, Graph) 
      { put(mis, v, MIS_UNFIX); }

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::maximal_independent_set<Graph, MISMap, IdDistribution,
						       PriorityQueueGenerator, MessageGenerator>
      D(g, mis, trans, idd, runtime_params.threads, sched_getcpu(), msg_gen);
    
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

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

#ifdef MIS_PRIORITY
#ifdef MIS_STATS
    D.print_pq_sizes();
#endif
#endif

#ifdef MIS_STATS
    D.print_stats();
#endif

    if (mis_params.verify) {
      if (trans.rank()==0)
	std::cout << "Verifying mis ..." << std::endl;

      if (!verify_mis(trans, g, mis)) {
	std::cout << "MIS Verification Failed" << std::endl;
	assert(false);
	return 0;
      }
    }

    vertices_size_type visited = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (get(mis, v) != MIS_UNFIX)
	++visited; 
    }
	  
    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
      r(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r(visited);
	  
    if (mis_params.verify)
      if (trans.rank() == 0)
	std::cout << "Visited " << total << " vertices of " << gparams.n << " in " << print_time(end - start) 
		  << std::endl;

    //if (total < 100) return -1.;
	  
    return end - start;
  }


  // FIX-Bucket algorithm
  template <typename Graph, typename IdDistribution,
	    typename graph_create_params,
	    typename MessageGenerator>
  time_type
  run_fix_mis_bucket(amplusplus::transport& trans, 
		     amplusplus::transport& barrier_trans, 
		     Graph& g,  
		     const IdDistribution& idd,
		     MessageGenerator msg_gen, 
		     graph_create_params& gparams,
		     instance_params& runtime_params,
		     mis_instance_params& mis_params) {

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;

    std::vector<state_t> misvec(num_vertices(g), MIS_UNFIX);
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    typedef boost::iterator_property_map<typename std::vector<state_t>::iterator, VertexIndexMap>  MISMap;
    MISMap mis(misvec.begin(), get(boost::vertex_index, g));

    if (trans.rank() == 0)
      std::cout << "Initializing mis-delta map ..." << std::endl;

    BGL_FORALL_VERTICES_T(v, g, Graph) 
      { put(mis, v, MIS_UNFIX); }

    trans.set_nthreads(runtime_params.threads);

    if (trans.rank() == 0)
      std::cout << "Creating algorithm instance ..." << std::endl;

    boost::graph::distributed::maximal_independent_set_delta<Graph, MISMap, IdDistribution,
							     append_buffer<Vertex, 10u>, MessageGenerator> 
      D(g, mis, trans, idd, runtime_params.flush, sched_getcpu(), msg_gen);
    
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

#ifdef MIS_STATS
    D.print_stats();
#endif

    // Back to one thread
    trans.set_nthreads(1);
    clear_thread_core_data();

    if (mis_params.verify) {
      if (trans.rank()==0)
	std::cout << "Verifying delta mis ..." << std::endl;

      if (!verify_mis(trans, g, mis)) {
	std::cout << "Bucket-MIS Verification Failed" << std::endl;
	assert(false);
	return 0;
      }
    }

    vertices_size_type visited = 0;
    BGL_FORALL_VERTICES_T(v, g, Graph) { 
      if (get(mis, v) != MIS_UNFIX)
	++visited; 
    }
	  
    boost::parallel::all_reduce<vertices_size_type, std::plus<vertices_size_type> > 
      r(trans, std::plus<vertices_size_type>());
    vertices_size_type total = r(visited);
	  
    if (mis_params.verify)
      if (trans.rank() == 0)
	std::cout << "Visited " << total << " vertices of " << gparams.n << " in " << print_time(end - start) 
		  << std::endl;

    return end - start;
  }


public:
  template <typename Graph, typename MessageGenerator,
	    typename graph_create_params>
  time_type operator()(const Graph& g, 
		       amplusplus::transport& trans, 
		       MessageGenerator& msg_gen,
		       graph_create_params& gparams,
		       instance_params& runtime_params,
		       mis_instance_params& mis_params) { // remove

    if (mis_params.algorithm == fix_bucket) {
      amplusplus::transport barrier_trans = trans.clone();
      if (mis_params.id_distribution == vertical) {
	return run_fix_mis_bucket(trans,
				  barrier_trans,
				  g,
				  block_id_distribution<Graph>(g, gparams.n),
				  msg_gen,
				  gparams,
				  runtime_params,
				  mis_params); 
      } else if (mis_params.id_distribution == horizontal) {
	return run_fix_mis_bucket(trans,
				  barrier_trans,
				  g,
				  row_id_distribution<Graph>(g, trans.size()),
				  msg_gen,
				  gparams,
				  runtime_params,
				  mis_params); 

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
				   

    }

    if (mis_params.algorithm == fix) {
      amplusplus::transport barrier_trans = trans.clone();

      if (mis_params.id_distribution == vertical) {
	return run_fix_mis(trans,
			   barrier_trans,
			   g,
			   block_id_distribution<Graph>(g, gparams.n),
			   msg_gen,
			   gparams,
			   runtime_params,
			   mis_params); //remove
      } else if (mis_params.id_distribution == horizontal) {
	return run_fix_mis(trans,
			   barrier_trans,
			   g,
			   row_id_distribution<Graph>(g, trans.size()),
			   msg_gen,
			   gparams,
			   runtime_params,
			   mis_params); //remove

      } else {
	std::cerr << "Invalid id distribution ! " << std::endl;
	assert(false);
      }
    }


    if (mis_params.algorithm == luby_a) {
      amplusplus::transport barrier_trans = trans.clone();
      return run_luby_maximal_is<boost::graph::distributed::select_a_functor_gen>(trans,
										  barrier_trans,
										  g,
										  msg_gen,
										  gparams,
										  runtime_params,
										  mis_params);

    }

    if (mis_params.algorithm == luby_av1) {
      amplusplus::transport barrier_trans = trans.clone();
      return run_luby_maximal_is<boost::graph::distributed::select_a_vertex_functor_gen>(trans,
										  barrier_trans,
										  g,
										  msg_gen,
										  gparams,
										  runtime_params,
										  mis_params);

    }


    if (mis_params.algorithm == luby_av2) {
      amplusplus::transport barrier_trans = trans.clone();
      return run_luby_maximal_is<boost::graph::distributed::select_a_v2_functor_gen>(trans,
										  barrier_trans,
										  g,
										  msg_gen,
										  gparams,
										  runtime_params,
										  mis_params);

    }

    if (mis_params.algorithm == luby_b) {
      amplusplus::transport barrier_trans = trans.clone();
      return run_luby_maximal_is<boost::graph::distributed::select_b_functor_gen>(trans,
										  barrier_trans,
										  g,
										  msg_gen,
										  gparams,
										  runtime_params,
										  mis_params);

    }




    
  }

};


int main(int argc, char* argv[]) {
  std::cout << "printing core id for process ..." << std::endl;
  print_core_id();

  executor<MISExecutor, mis_params, mis_instance_params> mis_executor;
  mis_executor.execute(argc, argv);  
  
}
