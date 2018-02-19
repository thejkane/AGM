// Copyright (C) 2004-2012 The Trustees of Indiana University.
// Pure Shiloach-Vishkin Connected Components
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Amila Kanewala
//           Nicholas Edmonds
//           Andrew Lumsdaine
//=========================================================
// Pure Shiloach-Vishkin algorithm without any optmizations
//=========================================================

#ifndef BOOST_GRAPH_PARALLEL_SV_CC_HPP
#define BOOST_GRAPH_PARALLEL_SV_CC_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

//#define BFS_SV_CC_HACK

#include <boost/property_map/parallel/lock_map.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/property_map/parallel/caching_property_map.hpp>
#include <boost/graph/parallel/algorithm.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/overloading.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/graph/parallel/properties.hpp>
#include <boost/graph/distributed/local_subgraph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/parallel/append_buffer.hpp>
#include <boost/optional.hpp>
#include <algorithm>
#include <vector>
#include <list>
#include <atomic>
#include <boost/graph/parallel/container_traits.hpp>
#include <boost/graph/parallel/iteration_macros.hpp>

#include <am++/counter_coalesced_message_type.hpp>

namespace boost { namespace graph { namespace distributed {

template <typename Graph, typename ParentMap,
            typename MessageGenerator = 
              amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen>,
            typename LockMap = boost::parallel::lock_map<typename property_map<Graph, vertex_index_t>::type>,
            typename VertexCompare = cc_vertex_compare<typename property_map<Graph, vertex_owner_t>::const_type,
                                                       typename property_map<Graph, vertex_local_t>::const_type>,
            typename VertexCombine = cc_vertex_combine<Graph,
                                                       typename property_map<Graph, vertex_owner_t>::const_type,
                                                       typename property_map<Graph, vertex_local_t>::const_type> >
class sv_cc {

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  // TODO: Add check that property_traits<ParentMap>::value_type is in fact Vertex

  typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap;

  typedef std::pair<Vertex, Vertex> vertex_pair_data;
  struct pointer_double_handler;

  struct hook_handler;
  struct update_parent_handler;

  // This message requests the recipient to attempt to hook p(u) onto p(v)
  typedef typename MessageGenerator::template call_result<vertex_pair_data, hook_handler, 
							  vertex_pair_owner<OwnerMap>,
							  amplusplus::idempotent_combination_t<VertexCombine, Vertex> >::type
  HookMessage;
  // duplicate removal is the best we can do here since we can't know p(p(v)

  // This message sends a request to the current parent to update a vertex's parent pointer
  // TODO: We could probably keep a flag indicating whether the parent's parent had changed 
  //       and only send back an UpdateParentMessage when needed
  typedef typename MessageGenerator::template call_result<vertex_pair_data, pointer_double_handler, 
							  vertex_pair_owner<OwnerMap>, amplusplus::no_reduction_t>::type
  PointerDoubleMessage;

  // This message is the response to PointerDoubleMessage above
  typedef typename MessageGenerator::template call_result<vertex_pair_data, update_parent_handler, 
							  vertex_pair_owner<OwnerMap>, amplusplus::no_reduction_t >::type
  UpdateParentMessage;

public:

  // TODO: ctor that takes a lock map

  sv_cc(amplusplus::transport& trans, Graph& g, const ParentMap& parent, LockMap& locks,
	MessageGenerator message_gen = 
	MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_pair_data>(), 0)),
      transport(trans), g(g), parent(parent), owner(get(vertex_owner, g)), 
      locks(locks), compare(get(vertex_owner, g), get(vertex_local, g)), 
      combine(get(vertex_owner, g), get(vertex_local, g)), level_sync(false), 
      hook_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner),
	       amplusplus::idempotent_combination(combine, graph_traits<Graph>::null_vertex())),
      pointer_double_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction),
      update_parent_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction)
  {
    initialize();
  }

  sv_cc(amplusplus::transport& trans, Graph& g, const ParentMap& parent, LockMap& locks, 
	VertexCompare compare, VertexCombine combine,
	MessageGenerator message_gen = 
	MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_pair_data>(), 0)),
      transport(trans), g(g), parent(parent), owner(get(vertex_owner, g)), 
      locks(locks), compare(compare), combine(combine), level_sync(false),
      hook_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner),
	       amplusplus::idempotent_combination(combine, graph_traits<Graph>::null_vertex())),
      pointer_double_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction),
      update_parent_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction)
  {
    initialize();
  }

  sv_cc(amplusplus::transport& trans, Graph& g, const ParentMap& parent,
	MessageGenerator message_gen = 
	MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_pair_data>(), 0)),
      transport(trans), g(g), parent(parent), owner(get(vertex_owner, g)), 
      locks(LockMap(get(vertex_index, g), 1 /*num_vertices(g)*/)), compare(get(vertex_owner, g), get(vertex_local, g)), 
      combine(get(vertex_owner, g), get(vertex_local, g)), level_sync(false),
      hook_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner),
	       amplusplus::idempotent_combination(combine, graph_traits<Graph>::null_vertex())),
      pointer_double_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction),
      update_parent_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction)
  {
    initialize();
  }

  sv_cc(amplusplus::transport& trans, Graph& g, const ParentMap& parent, VertexCompare compare, VertexCombine combine,
	MessageGenerator message_gen = 
	MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_pair_data>(), 0)),
      transport(trans), g(g), parent(parent), owner(get(vertex_owner, g)), 
      locks(LockMap(get(vertex_index, g), 1 /*num_vertices(g)*/)), compare(compare), combine(combine), level_sync(false),
      hook_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), 
	       amplusplus::idempotent_combination(combine, graph_traits<Graph>::null_vertex())),
      pointer_double_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction),
      update_parent_msg(message_gen, transport, vertex_pair_owner<OwnerMap>(owner), amplusplus::no_reduction)
  {
    initialize();
  }

  void set_level_sync() { // force level-synchronized hooking
    level_sync = true; 
  } 

  void initialize() {
    hook_msg.set_handler(hook_handler(*this));
    pointer_double_msg.set_handler(pointer_double_handler(*this));
    update_parent_msg.set_handler(update_parent_handler(*this));

#ifdef CC_PROFILING
    std::cout << "initializing called" << std::endl;
    pv_pointer_double_msgs.store(0);
    ppv_pointer_double_msgs.store(0);
    hook_msgs.store(0);
    local_hooks.store(0);
#endif
  }


  void print_stats() {
#ifdef CC_PROFILING
    uint64_t pv_pd = pv_pointer_double_msgs.load();
    uint64_t ppv_pd = ppv_pointer_double_msgs.load();
    uint64_t hk_msgs = hook_msgs.load();
    uint64_t lcl_hks = local_hooks.load();

    std::cout << "before all reduce pv_pd=" << pv_pd << std::endl;
    std::cout << "before all reduce ppv_pd=" << ppv_pd << std::endl;

    pv_pd = boost::parallel::all_reduce_sum(transport, pv_pd);
    ppv_pd = boost::parallel::all_reduce_sum(transport, ppv_pd);
    hk_msgs = boost::parallel::all_reduce_sum(transport, hk_msgs);
    lcl_hks = boost::parallel::all_reduce_sum(transport, lcl_hks);

    uint64_t pointer_dbls = pv_pd + ppv_pd;

    if (transport.rank() == 0) {
      std::cout << "Total pointer double messages : " << pointer_dbls
		<< ", first parent pointer doubles : " << pv_pd
		<< ", second parent pointer doubles : " << ppv_pd
		<< ", total hook messages : " << hk_msgs
		<< ", local hooks : " << lcl_hks
		<< std::endl;
    }
#endif

  }

  void operator() (unsigned int tid) { run(tid); }
  void run(unsigned int tid = 0);

  template <typename ComponentMap>
    typename property_traits<ComponentMap>::value_type
    number_components(const ComponentMap& component) {

    typedef typename property_traits<ComponentMap>::value_type component_size_type;

    typedef typename property_traits<ParentMap>::value_type parent_type;

    std::vector<parent_type> roots;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      parent_type pv = get(parent, v);

      // Null vertex is valid for parallel-search + SV
      // assert(pv != graph_traits<Graph>::null_vertex());

      if (std::find(roots.begin(), roots.end(), pv) == roots.end())
	roots.push_back(pv);
    }

    std::vector<parent_type> all_roots;

    boost::parallel::all_gather<parent_type> root_gather(transport);
    root_gather(roots, all_roots);

    std::sort(all_roots.begin(), all_roots.end(), compare);
    typename std::vector<parent_type>::iterator new_end = unique(all_roots.begin(), all_roots.end());

    std::map<parent_type, int> component_numbers;
    int i = 0;
    for (typename std::vector<parent_type>::iterator iter = all_roots.begin(); iter != new_end; ++iter, ++i)
      component_numbers[*iter] = i;

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      parent_type pv = get(parent, v);

      assert(component_numbers.find(pv) != component_numbers.end());

      put(component, v, component_numbers[pv]);
    }

    return static_cast<component_size_type>(new_end - all_roots.begin());
  }

private:

  bool try_hook(Vertex v, Vertex u) {
    Vertex pv;

    bool hooked;
    do {
      pv = get(parent, v);
      if (!compare(u, pv)) // if u is not better than v's current parent
	return false;
        
      // TODO: if (nthreads == 1) { put(parent, v, u); return true; }
      typename LockMap::value_type lock = locks.template maybe_lock<boost::parallel::atomics_supported<Vertex> >(v);
      hooked = maybe_exchange(parent, v, pv, u);
    } while(!hooked);

    return true;
  }

  /*  template <typename LocalIndexMap>  // property_traits<LocalIndexMap>::key_type is Vertex, but we can only use template parameters of the function in enable_if
  typename boost::enable_if<boost::parallel::atomics_supported<typename property_traits<LocalIndexMap>::key_type>, void>::type
  find_local_roots(LocalIndexMap& local_index, int tid, int nthreads) {
    BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {

#ifdef BFS_SV_CC_HACK
      if (get(parent, v) == graph_traits<Graph>::null_vertex()) continue;
#endif
      vertices_size_type component = local_subgraph_component[get(local_index, v)];
      Vertex current_root;
        
      using boost::parallel::bool_compare_and_swap;
        
      do {

	current_root = roots[component];
	if (current_root != graph_traits<Graph>::null_vertex() &&
	    get(local_index, v) >= get(local_index, current_root))
	  break;
      } while(!bool_compare_and_swap(&roots[component], current_root, v));
      // If this vertex is smaller than the current root of it's component
      // try to make it the root
    }
    }*/

  /*template <typename LocalIndexMap>
  typename boost::disable_if<boost::parallel::atomics_supported<typename property_traits<LocalIndexMap>::key_type> >::type
  find_local_roots(LocalIndexMap& local_index, int tid, int) {

    if (tid == 0)
      BGL_FORALL_VERTICES_T(v, g, Graph) {
#ifdef BFS_SV_CC_HACK
	if (get(parent, v) == graph_traits<Graph>::null_vertex()) continue;
#endif
	vertices_size_type component = local_subgraph_component[get(local_index, v)];
          
	if (roots[component] == graph_traits<Graph>::null_vertex() ||
	    get(local_index, v) < get(local_index, roots[component]))
	  roots[component] = v;
      }
      }*/

  const int dummy_first_member_for_init_order; // Unused

  amplusplus::transport& transport;
  Graph&                 g;
  const ParentMap&       parent;
  const OwnerMap         owner;
  LockMap&               locks;
  VertexCompare          compare;
  VertexCombine          combine; 
  bool                   level_sync;

  // For level-synchronized version
  append_buffer<std::pair<Vertex, Vertex> > hook_requests;

  shared_ptr<amplusplus::detail::barrier> t_bar;

  HookMessage          hook_msg;
  PointerDoubleMessage pointer_double_msg;
  UpdateParentMessage  update_parent_msg;

  // Shared data structures for threading
  std::vector<Vertex> roots;
  std::vector<vertices_size_type> local_subgraph_component; // TODO: Perhaps rather than this vector we should store a property map of the local subgraph

  long unsigned int hooked, doubled;

  boost::scoped_ptr<local_subgraph<Graph> > ls;
#ifdef CC_PROFILING
  std::atomic<std::uint64_t> pv_pointer_double_msgs;
  std::atomic<std::uint64_t> ppv_pointer_double_msgs;
  std::atomic<std::uint64_t> hook_msgs;
  std::atomic<std::uint64_t> local_hooks;
#endif
};

#define SV_CC_PARMS                                   \
  typename Graph, typename ParentMap, typename MessageGenerator, typename LockMap, typename VertexCompare, typename VertexCombine

#define SV_CC_TYPE                                    \
  sv_cc<Graph, ParentMap, MessageGenerator, LockMap, VertexCompare, VertexCombine>

template <SV_CC_PARMS>
void
SV_CC_TYPE::run(unsigned int tid)
{
  AMPLUSPLUS_WITH_THREAD_ID(tid) {

    amplusplus::transport::rank_type rank = transport.rank();
    //std::cout << "Starting thread - " << tid << " on rank : "
    //	      << rank << std::endl;
    size_t nthreads = transport.get_nthreads();
    
    if (tid == 0)
      t_bar.reset(new amplusplus::detail::barrier(nthreads));

    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized 
    { amplusplus::scoped_epoch epoch(transport); } 

#ifdef PBGL_TIMING
    double start = MPI_Wtime();
#endif

    // Parallel phase
    //
    do {
      t_bar->wait(); // Make sure we don't change hooked before all threads use it in the while() below
      if (tid == 0) hooked = 0;
      t_bar->wait();
      
      {
        amplusplus::scoped_epoch epoch(transport);
        
	BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
          Vertex pv = get(parent, v);
          
#ifdef BFS_SV_CC_HACK
          if (pv == graph_traits<Graph>::null_vertex()) continue;  // This is a hack for BFS + SV CC 
#endif
	  std::set<Vertex> parallel_edges; // to avoid parallel edges
	  BGL_FORALL_ADJ_T(v, u, g, Graph) {
#ifdef CC_PROFILING
  	    if (get(owner, u) != rank) {
	      hook_msgs++;
	    } else
	      local_hooks++;
#endif	 
	    if (v != u) { // to avoid self-loops
	      if (parallel_edges.insert(u).second)
		hook_msg.send(vertex_pair_data(u, pv));
	    }
 	  }
	}
      }
      
      if (level_sync) {
        amplusplus::scoped_epoch epoch(transport);

	// Process hook_msgs
	for (int i = tid ; i < hook_requests.size() ; i += nthreads)
	  hook_msg.send(hook_requests[i]);

        t_bar->wait();
        if (tid == 0) hook_requests.clear();
        t_bar->wait();
      }

      t_bar->wait(); // alternately call wait_for_epoch_cleanup() above
      // Can't use scoped_epoch_value because value may be written by handlers in TD
      hooked = boost::parallel::all_reduce_sum(transport, hooked);
      
#ifdef PRINT_DEBUG
      if (tid == 0 && transport.rank() == 0) std::cout << "Hooked: " << hooked << std::endl; 
#endif
      
      do {
        t_bar->wait();
        if (tid == 0) doubled = 0;
        t_bar->wait();
        
        {
          amplusplus::scoped_epoch epoch(transport);
          
	  BGL_PARFORALL_VERTICES_T(v, g, Graph, tid, nthreads) {
	    //BGL_FORALL_VERTICES_T(v, g, Graph) {
            Vertex pv = get(parent, v);
            
#ifdef BFS_SV_CC_HACK
            if (pv == graph_traits<Graph>::null_vertex()) continue; // This is a hack for BFS + SV CC 
#endif
	    if (get(owner, pv) == rank) {
	      Vertex ppv = get(parent, pv);
	      if (ppv != pv) {
		doubled = 1;
		put(parent, v, ppv);
	      }
	    } else { // p(v) is remote
#ifdef CC_PROFILING
	      pv_pointer_double_msgs++;
#endif	    
	      pointer_double_msg.send(vertex_pair_data(pv, v));
	    }
	  }
	}
	
	t_bar->wait(); // alternately call wait_for_epoch_cleanup() above
	// Can't use scoped_epoch_value because value may be written by handlers in TD
	doubled = boost::parallel::all_reduce_sum(transport, doubled);
	
#ifdef PRINT_DEBUG
        if (tid == 0 && transport.rank() == 0) std::cout << "Doubled: " << doubled << std::endl; 
#endif
        
      } while (doubled > 0);

    } while (hooked > 0);

  } // End AMPLUSPLUS_WITH_THREAD_ID

  // HOOKING
  //
  // Find cluster root and hook it to minimum adjacent vertex of all vertices in component
  //
  // WHILE(some vertex hooks successfully)
  //
  // BGL_FORALL_VERTICES_T(v, g, Graph)
  //   BGL_FORALL_ADJ_T(v, u, g, Graph)
  //     if (u is remote)
  //       try_hook(u, p(v) // try to hook p(u) to p(v)
  //
  // Recipient of try_hook -- either p(u) is local and we attempt to hook or p(u) 
  // is remote and the message gets forwarded... this is better than the case where
  // all hooking requests take two messages as long as there is at least one local root

  // Owner of u responts with hook(p(v), p(u))
  //   which hooks the component v is in to the one u is in... there will be many of these 
  //   but only the lowest will succeed
  //   (Reduction op reduces write to same vertex (remote component root))

  // POINTER DOUBLING
  //
  // do {
  //   BGL_FORALL_VERTICES_T(v, g, Graph)
  //     p(v) = p(p(v));
  // } while (some p(v) changed);
}

template <SV_CC_PARMS>
struct SV_CC_TYPE::
hook_handler {

  hook_handler() : self(NULL) {}
  hook_handler(sv_cc& self) : self(&self) {}

  // Try to hook p(data.first) onto data.second
  void operator() (const vertex_pair_data& data) const
  {
    Vertex v = data.second;
    Vertex u = data.first;
    Vertex pu = get(self->parent, u);

    assert(v != graph_traits<Graph>::null_vertex());
    assert(u != graph_traits<Graph>::null_vertex());
    assert(pu != graph_traits<Graph>::null_vertex());

    // It would be possible to avoid this ownership check with a second
    // message type that attempts to hook unconditionally, but that would
    // require a whole additional message type...
    if (get(self->owner, pu) == self->transport.rank()) {
      if (self->try_hook(pu, v))
        self->hooked = 1; // Unconditional write, no synchronization required 
                          // and we don't have to worry about overflow on BG/P
    } else {
#ifdef CC_PROFILING
      self->hook_msgs++;
#endif
      if (self->level_sync) 
        self->hook_requests.push_back(std::make_pair(pu, v));
      else
        self->hook_msg.send(vertex_pair_data(pu, v)); // pu is a root so p(pu) = pu
    }
  }

protected:
  sv_cc* self;
};

template <SV_CC_PARMS>
struct SV_CC_TYPE::
pointer_double_handler {

  pointer_double_handler() : self(NULL) {}
  pointer_double_handler(sv_cc& self) : self(&self) {}

  // (vertex requesting parent for, vertex to update)
  void operator() (const vertex_pair_data& data) const
  {
    assert(data.first != graph_traits<Graph>::null_vertex());
    assert(data.second != graph_traits<Graph>::null_vertex());

#ifdef CC_PROFILING
    amplusplus::transport::rank_type rank = self->transport.rank();
    if (get(self->owner, data.second) != rank) {
      self->ppv_pointer_double_msgs++;
    }
#endif

    self->update_parent_msg.send(vertex_pair_data(data.second, get(self->parent, data.first)));
  }

protected:
  sv_cc* self;
};

template <SV_CC_PARMS>
struct SV_CC_TYPE::
update_parent_handler {

  update_parent_handler() : self(NULL) {}
  update_parent_handler(sv_cc& self) : self(&self) {}

  // (vertex to update, new parent value)
  void operator() (const vertex_pair_data& data) const
  {
    Vertex v = data.first;
    Vertex pv = get(self->parent, v);
    Vertex ppv = data.second;

    assert(v != graph_traits<Graph>::null_vertex());
    assert(pv != graph_traits<Graph>::null_vertex());
    assert(ppv != graph_traits<Graph>::null_vertex());

    if (pv != ppv) {
      self->doubled = 1;
      put(self->parent, v, ppv); // each vertex only sends one pointer-doubling message so no need for atomic update
    }
  }

protected:
  sv_cc* self;
};

} } } // end namespace boost::graph::distributed

#endif // BOOST_GRAPH_PARALLEL_CC_HPP
