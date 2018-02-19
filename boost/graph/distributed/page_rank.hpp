// Copyright (C) 2004-2012 The Trustees of Indiana University.
// Copyright (C) 2002 Brad King and Douglas Gregor

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine
#ifndef BOOST_PARALLEL_GRAPH_PAGE_RANK_HPP
#define BOOST_PARALLEL_GRAPH_PAGE_RANK_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/parallel/iteration_macros.hpp>
#include <boost/graph/page_rank.hpp>

namespace boost { namespace graph { namespace distributed {

namespace detail {
  template<typename T>
  struct rank_accumulate_reducer {
    BOOST_STATIC_CONSTANT(bool, non_default_resolver = true);

    template<typename K>
    T operator()(const K&) const { return T(0); }

    template<typename K>
    T operator()(const K&, const T& x, const T& y) const { return x + y; }
  };

  template<typename Graph, typename RankMap, typename RankMap2, typename OwnerMap>
  void page_rank_step(const Graph& g, RankMap from_rank, RankMap2 to_rank, OwnerMap owner,
                      typename property_traits<RankMap>::value_type damping,
                      int tid, int nthreads)
  {
    typedef typename property_traits<RankMap>::value_type rank_type;

    // Set new rank maps 
    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
      put(to_rank, v, rank_type(1 - damping));
    }

    {
//       amplusplus::scoped_epoch epoch(g.transport());

      BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
	if (out_degree(u, g) != 0) {
	  rank_type u_rank_out = damping * get(from_rank, u) / out_degree(u, g);
	  BGL_FORALL_ADJ_T(u, v, g, Graph)
	    if (get(owner, v) == g.transport().rank())
	      put(to_rank, v, get(to_rank, v) + u_rank_out);
	      // NGE: Calling get(to_rank, v) in the non-local case will work,
              //      but it will create lots of unnecessary ghost cells
	    else
	      put(to_rank, v, u_rank_out); // TODO: atomic_put
	}
      }
    }

    // NGE: If we have multiple threads we can use cm_forward in DPMap put()
    //      instead and skip the flush here... 
    synchronize(to_rank); // Flush ghost cells
    synchronize(from_rank); // Reset ghost cells
  }
} // end namespace detail

template <typename Graph, typename RankMap, typename Done = n_iterations,
 	  typename MessageGenerator = 
	    amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class page_rank {

  typedef typename property_traits<RankMap>::value_type rank_type;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

public:

  // TODO: Make DPmaps copy constructible so we can create rank2 from rank
  page_rank(Graph& g, RankMap rank, RankMap rank2, Done done, rank_type damping, 
	    vertices_size_type n)
    : g(g), transport(g.transport()), rank(rank), rank2(rank2), done(done), 
      damping(damping), n(n)
  {
    rank.set_reduce(detail::rank_accumulate_reducer<rank_type>());
    rank2.set_reduce(detail::rank_accumulate_reducer<rank_type>());
    rank.set_consistency_model(boost::parallel::cm_flush | boost::parallel::cm_reset);
    rank2.set_consistency_model(boost::parallel::cm_flush | boost::parallel::cm_reset);
  }

  void operator() (int tid = 0) {

    AMPLUSPLUS_WITH_THREAD_ID(tid) {

      size_t nthreads = transport.get_nthreads();

      if (tid == 0)
	t_bar.reset(new amplusplus::detail::barrier(nthreads));
      
      // This barrier acts as a temporary barrier until we can be sure t_bar is initialized 
      { amplusplus::scoped_epoch epoch(transport); } 

      rank_type initial_rank = rank_type(rank_type(1) / n);
      BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) put(rank, v, initial_rank);

      typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap;
      OwnerMap owner = get(vertex_owner, g);

      bool to_map_2 = true;
      while ((to_map_2 && !done(rank, g)) || (!to_map_2 && !done(rank2, g))) {

	/**
	 * PageRank can implemented slightly more efficiently on a
	 * bidirectional graph than on an incidence graph. However,
	 * distributed PageRank requires that we have the rank of the
	 * source vertex available locally, so we force the incidence
	 * graph implementation, which pushes rank from source to
	 * target.
	 */
	if (to_map_2) 
	  detail::page_rank_step(g, rank, rank2, owner, damping, tid, nthreads);
	else 
	  detail::page_rank_step(g, rank2, rank, owner, damping, tid, nthreads);

	to_map_2 = !to_map_2;
      }
      
      rank.reset();
      
      if (!to_map_2)
	BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) put(rank, v, get(rank2, v));
    }
  }

private:

  const Graph& g;
  amplusplus::transport& transport;
  RankMap rank, rank2;
  Done done;	    
  rank_type damping;
  vertices_size_type n;

  shared_ptr<amplusplus::detail::barrier> t_bar;
};

} } } // end namespace boost::graph::distributed

#endif // BOOST_PARALLEL_GRAPH_PAGE_RANK_HPP
