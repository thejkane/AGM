// Copyright (C) 2004-2012 The Trustees of Indiana University.
// Copyright (C) 2007 Douglas Gregor 

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Nicholas Edmonds
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_DISTRIBUTED_ADJACENCY_LIST_HPP
#define BOOST_GRAPH_DISTRIBUTED_ADJACENCY_LIST_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

// AM++
#include "am++/basic_coalesced_message_type.hpp"
#include <am++/make_mpi_datatype.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/parallel/local_property_map.hpp>
#include <boost/graph/parallel/detail/property_holders.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <cassert>
#include <list>
#include <algorithm>
#include <boost/limits.hpp>
#include <boost/graph/parallel/properties.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/distributed/selector.hpp>

#include <boost/graph/distributed/shuffled_distribution.hpp>

namespace boost {

  /// The type used to store an identifier that uniquely names a processor.
  typedef amplusplus::transport::rank_type rank_type;

  // Tell which processor the target of an edge resides on (for
  // directed graphs) or which processor the other end point of the
  // edge resides on (for undirected graphs).
  enum edge_target_processor_id_t { edge_target_processor_id };
  BOOST_INSTALL_PROPERTY(edge, target_processor_id);

  // For undirected graphs, tells whether the edge is locally owned.
  enum edge_locally_owned_t { edge_locally_owned };
  BOOST_INSTALL_PROPERTY(edge, locally_owned);

  // For bidirectional graphs, stores the incoming edges.
  enum vertex_in_edges_t { vertex_in_edges };
  BOOST_INSTALL_PROPERTY(vertex, in_edges);

  /// Tag class for directed, distributed adjacency list
  struct directed_distributed_adj_list_tag
    : public virtual distributed_graph_tag,
      public virtual distributed_vertex_list_graph_tag,
      public virtual distributed_edge_list_graph_tag,
      public virtual incidence_graph_tag,
      public virtual adjacency_graph_tag {};

  /// Tag class for bidirectional, distributed adjacency list
  struct bidirectional_distributed_adj_list_tag
    : public virtual distributed_graph_tag,
      public virtual distributed_vertex_list_graph_tag,
      public virtual distributed_edge_list_graph_tag,
      public virtual incidence_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual bidirectional_graph_tag {};

  /// Tag class for undirected, distributed adjacency list
  struct undirected_distributed_adj_list_tag
    : public virtual distributed_graph_tag,
      public virtual distributed_vertex_list_graph_tag,
      public virtual distributed_edge_list_graph_tag,
      public virtual incidence_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual bidirectional_graph_tag {};

  namespace detail { namespace parallel {
  
    /**
     * A distributed vertex descriptor. These descriptors contain both
     * the ID of the processor that owns the vertex and a local vertex
     * descriptor that identifies the particular vertex for that
     * processor.
     */
    template<typename LocalDescriptor>
    struct global_descriptor
    {
      typedef LocalDescriptor local_descriptor_type;

      global_descriptor() : owner(), local() { }

      // NGE: This is a hack for AM++ reductions
      //      Local descriptor should really be null_vertex() but I don't know how to get one here
      global_descriptor(int) : owner(), local(std::numeric_limits<LocalDescriptor>::max()) {}

      global_descriptor(rank_type owner, LocalDescriptor local)
        : owner(owner), local(local) { }

      rank_type owner;
      LocalDescriptor local;

      /**
       * A function object that, given a processor ID, generates
       * distributed vertex descriptors from local vertex
       * descriptors. This function object is used by the
       * vertex_iterator of the distributed adjacency list.
       */
      struct generator
      {
        typedef global_descriptor<LocalDescriptor> result_type;
        typedef LocalDescriptor argument_type;

        generator() {}
        generator(rank_type owner) : owner(owner) {}

        result_type operator()(argument_type v) const
        { return result_type(owner, v); }

      private:
        rank_type owner;
      };
    };

    /// Determine the process that owns the given descriptor
    template<typename LocalDescriptor>
    inline rank_type owner(const global_descriptor<LocalDescriptor>& v)
    { return v.owner; }

    /// Determine the local portion of the given descriptor
    template<typename LocalDescriptor>
    inline LocalDescriptor local(const global_descriptor<LocalDescriptor>& v)
    { return v.local; }

    /// Compare distributed vertex descriptors for equality
    template<typename LocalDescriptor>
    inline bool
    operator==(const global_descriptor<LocalDescriptor>& u,
               const global_descriptor<LocalDescriptor>& v)
    {
      return u.owner == v.owner && u.local == v.local;
    }

    /// Compare distributed vertex descriptors for inequality
    template<typename LocalDescriptor>
    inline bool
    operator!=(const global_descriptor<LocalDescriptor>& u,
               const global_descriptor<LocalDescriptor>& v)
    { return !(u == v); }

    template<typename LocalDescriptor>
    inline bool
    operator<(const global_descriptor<LocalDescriptor>& u,
              const global_descriptor<LocalDescriptor>& v)
    {
      return (u.owner) < v.owner || (u.owner == v.owner && (u.local) < v.local);
    }

    template<typename LocalDescriptor>
    inline bool
    operator<=(const global_descriptor<LocalDescriptor>& u,
               const global_descriptor<LocalDescriptor>& v)
    {
      return u.owner <= v.owner || (u.owner == v.owner && u.local <= v.local);
    }

    template<typename LocalDescriptor>
    inline bool
    operator>(const global_descriptor<LocalDescriptor>& u,
              const global_descriptor<LocalDescriptor>& v)
    {
      return v < u;
    }

    template<typename LocalDescriptor>
    inline bool
    operator>=(const global_descriptor<LocalDescriptor>& u,
               const global_descriptor<LocalDescriptor>& v)
    {
      return v <= u;
    }

    // DPG TBD: Add <, <=, >=, > for global descriptors

    /**
     * A Readable Property Map that extracts a global descriptor pair
     * from a global_descriptor.
     */
    template<typename LocalDescriptor>
    struct global_descriptor_property_map
    {
      typedef std::pair<rank_type, LocalDescriptor> value_type;
      typedef value_type reference;
      typedef global_descriptor<LocalDescriptor> key_type;
      typedef readable_property_map_tag category;
    };

    template<typename LocalDescriptor>
    inline std::pair<rank_type, LocalDescriptor>
    get(global_descriptor_property_map<LocalDescriptor>,
        global_descriptor<LocalDescriptor> x)
    {
      return std::pair<rank_type, LocalDescriptor>(x.owner, x.local);
    }

    /**
     * A Readable Property Map that extracts the owner of a global
     * descriptor.
     */
    template<typename LocalDescriptor>
    struct owner_property_map
    {
      typedef rank_type value_type;
      typedef value_type reference;
      typedef global_descriptor<LocalDescriptor> key_type;
      typedef readable_property_map_tag category;
    };

    template<typename LocalDescriptor>
    inline rank_type
    get(owner_property_map<LocalDescriptor>,
        global_descriptor<LocalDescriptor> x)
    {
      return x.owner;
    }

    /**
     * A Readable Property Map that extracts the local descriptor from
     * a global descriptor.
     */
    template<typename LocalDescriptor>
    struct local_descriptor_property_map
    {
      typedef LocalDescriptor value_type;
      typedef value_type reference;
      typedef global_descriptor<LocalDescriptor> key_type;
      typedef readable_property_map_tag category;
    };

    template<typename LocalDescriptor>
    inline LocalDescriptor
    get(local_descriptor_property_map<LocalDescriptor>,
        global_descriptor<LocalDescriptor> x)
    {
      return x.local;
    }

    /**
     * Stores an incoming edge for a bidirectional distributed
     * adjacency list. The user does not see this type directly,
     * because it is just an implementation detail.
     */
    template<typename Edge>
    struct stored_in_edge
    {
      stored_in_edge(rank_type sp, Edge e)
        : source_processor(sp), e(e) {}

      rank_type source_processor;
      Edge e;
    };

    /**
     * A distributed edge descriptor. These descriptors contain the
     * underlying edge descriptor, the processor IDs for both the
     * source and the target of the edge, and a boolean flag that
     * indicates which of the processors actually owns the edge.
     */
    template<typename Edge>
    struct edge_descriptor
    {
      edge_descriptor(rank_type sp = rank_type(),
                      rank_type tp = rank_type(),
                      Edge ld = Edge())
        : source_processor(sp), target_processor(tp), local(ld) {}

      rank_type owner() const
      {
        /* return source_owns_edge? source_processor : target_processor; */
	return source_processor;
      }

      /// The processor that the source vertex resides on
      rank_type source_processor;

      /// The processor that the target vertex resides on
      rank_type target_processor;

      /// The local edge descriptor.
      Edge local;

      /**
       * Function object that generates edge descriptors for the
       * out_edge_iterator of the given distributed adjacency list
       * from the edge descriptors of the underlying adjacency list.
       */
      template<typename Graph>
      class out_generator
      {
        typedef typename Graph::directed_selector directed_selector;

      public:
        typedef edge_descriptor<Edge> result_type;
        typedef Edge argument_type;

        out_generator() : g(0) {}
        explicit out_generator(const Graph& g) : g(&g) {}

        result_type operator()(argument_type e) const
        { return map(e, directed_selector()); }

      private:
        result_type map(argument_type e, directedS) const
        {
          return result_type(g->processor(),
                             get(edge_target_processor_id, g->base(), e),
                             e);
        }

        const Graph* g;
      };

      /**
       * Function object that generates edge descriptors for the
       * in_edge_iterator of the given distributed adjacency list
       * from the edge descriptors of the underlying adjacency list.
       */
      template<typename Graph>
      class in_generator
      {
        typedef typename Graph::directed_selector DirectedS;

      public:
        typedef Edge argument_type;
        typedef edge_descriptor<Edge> result_type;

        in_generator() : g(0) {}
        explicit in_generator(const Graph& g) : g(&g) {}

        result_type operator()(argument_type e) const
        { return map(e, DirectedS()); }

      private:
        const Graph* g;
      };
    };

    /// Determine the process that owns this edge
    template<typename Edge>
    inline rank_type
    owner(const edge_descriptor<Edge>& e)
    /* { return e.source_owns_edge? e.source_processor : e.target_processor; } */
    {return e.source_processor; }

    /// Determine the local descriptor for this edge.
    template<typename Edge>
    inline Edge
    local(const edge_descriptor<Edge>& e)
    { return e.local; }

    /**
     * A Readable Property Map that extracts the owner and local
     * descriptor of an edge descriptor.
     */
    template<typename Edge>
    struct edge_global_property_map
    {
      typedef std::pair<rank_type, Edge> value_type;
      typedef value_type reference;
      typedef edge_descriptor<Edge> key_type;
      typedef readable_property_map_tag category;
    };

    template<typename Edge>
    inline std::pair<rank_type, Edge>
    get(edge_global_property_map<Edge>, const edge_descriptor<Edge>& e)
    {
      typedef std::pair<rank_type, Edge> result_type;
//       return result_type(e.source_owns_edge? e.source_processor
//                          /* target owns edge*/: e.target_processor,
//                          e.local);
      return result_type(e.source_processor, e.local);
    }

    /**
     * A Readable Property Map that extracts the owner of an edge
     * descriptor.
     */
    template<typename Edge>
    struct edge_owner_property_map
    {
      typedef rank_type value_type;
      typedef value_type reference;
      typedef edge_descriptor<Edge> key_type;
      typedef readable_property_map_tag category;
    };

    template<typename Edge>
    inline rank_type
    get(edge_owner_property_map<Edge>, const edge_descriptor<Edge>& e)
    /* { return e.source_owns_edge? e.source_processor : e.target_processor; } */
    { return e.source_processor; }

    /**
     * A Readable Property Map that extracts the local descriptor from
     * an edge descriptor.
     */
    template<typename Edge>
    struct edge_local_property_map
    {
      typedef Edge value_type;
      typedef value_type reference;
      typedef edge_descriptor<Edge> key_type;
      typedef readable_property_map_tag category;
    };

    template<typename Edge>
    inline Edge
    get(edge_local_property_map<Edge>,
        const edge_descriptor<Edge>& e)
    {
      return e.local;
    }

    /** Compare distributed edge descriptors for equality.
     *
     * \todo need edge_descriptor to know if it is undirected so we
     * can compare both ways.
     */
    template<typename Edge>
    inline bool
    operator==(const edge_descriptor<Edge>& e1,
               const edge_descriptor<Edge>& e2)
    {
      return (e1.source_processor == e2.source_processor
              && e1.target_processor == e2.target_processor
              && e1.local == e2.local);
    }

    /// Compare distributed edge descriptors for inequality.
    template<typename Edge>
    inline bool
    operator!=(const edge_descriptor<Edge>& e1,
               const edge_descriptor<Edge>& e2)
    { return !(e1 == e2); }

    /**
     * Configuration for the distributed adjacency list. We use this
     * parameter to store all of the configuration details for the
     * implementation of the distributed adjacency list, which allows us to
     * get at the distribution type in the maybe_named_graph.
     */
    template<typename OutEdgeListS, typename InVertexListS, 
	     typename InDistribution, typename DirectedS, 
	     typename VertexProperty, typename EdgeProperty, 
	     typename GraphProperty, typename EdgeListS>
    struct adjacency_list_config
    {
      typedef typename mpl::if_<is_same<InVertexListS, defaultS>, 
                                vecS, InVertexListS>::type 
        VertexListS;

      /// Introduce the target processor ID property for each edge
      typedef property<edge_target_processor_id_t, rank_type,
                       EdgeProperty> edge_property_with_id;

      /// For undirected graphs, introduce the locally-owned property for edges
      typedef edge_property_with_id base_edge_property_type;

      /// The edge descriptor type for the local subgraph
      typedef typename adjacency_list_traits<OutEdgeListS,
                                             VertexListS,
                                             directedS>::edge_descriptor
        local_edge_descriptor;

      // Bidirectional graphs have an extra vertex property to store
      // the incoming edges.
      typedef VertexProperty base_vertex_property_type;

      // The type of the distributed adjacency list
      typedef adjacency_list<OutEdgeListS,
                             distributedS<VertexListS, 
                                          InDistribution>,
                             DirectedS, VertexProperty, EdgeProperty,
                             GraphProperty, EdgeListS> 
        graph_type;

      // The type of the underlying adjacency list implementation
      typedef adjacency_list<OutEdgeListS, VertexListS, directedS,
                             base_vertex_property_type,
                             base_edge_property_type,
                             GraphProperty,
                             EdgeListS> 
        inherited;
      
      typedef InDistribution in_distribution_type;
      typedef typename inherited::vertices_size_type vertices_size_type;

        typedef typename mpl::if_<is_same<InDistribution, defaultS>,
				  boost::parallel::variant_distribution<vertices_size_type>,
				  InDistribution>::type base_distribution_type;

          typedef ::boost::graph::distributed::shuffled_distribution<
          base_distribution_type> distribution_type;

      typedef VertexProperty vertex_property_type;
      typedef EdgeProperty edge_property_type;

      typedef VertexListS vertex_list_selector;
      typedef OutEdgeListS out_edge_list_selector;
      typedef DirectedS directed_selector;
      typedef GraphProperty graph_property_type;
      typedef EdgeListS edge_list_selector;
    };

    // Maybe initialize the indices of each vertex
    template<typename IteratorPair, typename VertexIndexMap>
    void
    maybe_initialize_vertex_indices(IteratorPair p, VertexIndexMap to_index,
                                    read_write_property_map_tag)
    {
      typedef typename property_traits<VertexIndexMap>::value_type index_t;
      index_t next_index = 0;
      while (p.first != p.second)
        put(to_index, *p.first++, next_index++);
    }

    template<typename IteratorPair, typename VertexIndexMap>
    inline void
    maybe_initialize_vertex_indices(IteratorPair , VertexIndexMap ,
                                    readable_property_map_tag)
    {
      // Do nothing
    }

    template<typename IteratorPair, typename VertexIndexMap>
    inline void
    maybe_initialize_vertex_indices(IteratorPair p, VertexIndexMap to_index)
    {
      typedef typename property_traits<VertexIndexMap>::category category;
      maybe_initialize_vertex_indices(p, to_index, category());
    }

    template<typename IteratorPair>
    inline void
    maybe_initialize_vertex_indices(IteratorPair p,
                                    ::boost::detail::error_property_not_found)
    { }

    //------------------------------------------------------------------------
    // Distributed adjacency list property map details
    /**
     * Metafunction that extracts the given property from the given
     * distributed adjacency list type. This could be implemented much
     * more cleanly, but even newer versions of GCC (e.g., 3.2.3)
     * cannot properly handle partial specializations involving
     * enumerator types.
     */
    template<typename Property>
    struct get_adj_list_pmap
    {
      template<typename Graph>
      struct apply
      {
        typedef Graph graph_type;
        typedef typename graph_type::inherited base_graph_type;
        typedef typename property_map<base_graph_type, Property>::type
          local_pmap;
        typedef typename property_map<base_graph_type, Property>::const_type
          local_const_pmap;

        typedef graph_traits<graph_type> traits;
        typedef typename graph_type::local_vertex_descriptor local_vertex;
        typedef typename property_traits<local_pmap>::key_type local_key_type;

        typedef typename property_traits<local_pmap>::value_type value_type;

        typedef typename property_map<Graph, vertex_global_t>::const_type
          vertex_global_map;
        typedef typename property_map<Graph, edge_global_t>::const_type
          edge_global_map;

        typedef typename mpl::if_c<(is_same<local_key_type,
                                            local_vertex>::value),
                                   vertex_global_map, edge_global_map>::type
          global_map;

      public:
        typedef ::boost::parallel::distributed_property_map<
                  global_map, local_pmap> type;

        typedef ::boost::parallel::distributed_property_map<
                  global_map, local_const_pmap> const_type;
      };
    };

    /**
     * The local vertex index property map is actually a mapping from
     * the local vertex descriptors to vertex indices.
     */
    template<>
    struct get_adj_list_pmap<vertex_local_index_t>
    {
      template<typename Graph>
      struct apply
        : ::boost::property_map<typename Graph::inherited, vertex_index_t>
      { };
    };

    /**
     * The vertex index property map maps from global descriptors
     * (e.g., the vertex descriptor of a distributed adjacency list)
     * to the underlying local index. It is not valid to use this
     * property map with nonlocal descriptors.
     */
    template<>
    struct get_adj_list_pmap<vertex_index_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename property_map<Graph, vertex_global_t>::const_type
          global_map;

        typedef property_map<typename Graph::inherited, vertex_index_t> local;

      public:
        typedef local_property_map<global_map,
                                   typename local::type> type;
        typedef local_property_map<global_map,
                                   typename local::const_type> const_type;
      };
    };

    /**
     * The vertex owner property map maps from vertex descriptors to
     * the processor that owns the vertex.
     */
    template<>
    struct get_adj_list_pmap<vertex_global_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename Graph::local_vertex_descriptor
          local_vertex_descriptor;
      public:
        typedef global_descriptor_property_map<local_vertex_descriptor> type;
        typedef type const_type;
      };
    };

    /**
     * The vertex owner property map maps from vertex descriptors to
     * the processor that owns the vertex.
     */
    template<>
    struct get_adj_list_pmap<vertex_owner_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename Graph::local_vertex_descriptor
          local_vertex_descriptor;
      public:
        typedef owner_property_map<local_vertex_descriptor> type;
        typedef type const_type;
      };
    };

    /**
     * The vertex local property map maps from vertex descriptors to
     * the local descriptor for that vertex.
     */
    template<>
    struct get_adj_list_pmap<vertex_local_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename Graph::local_vertex_descriptor
          local_vertex_descriptor;
      public:
        typedef local_descriptor_property_map<local_vertex_descriptor> type;
        typedef type const_type;
      };
    };

    /**
     * The edge global property map maps from edge descriptors to
     * a pair of the owning processor and local descriptor.
     */
    template<>
    struct get_adj_list_pmap<edge_global_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename Graph::local_edge_descriptor
          local_edge_descriptor;
      public:
        typedef edge_global_property_map<local_edge_descriptor> type;
        typedef type const_type;
      };
    };

    /**
     * The edge owner property map maps from edge descriptors to
     * the processor that owns the edge.
     */
    template<>
    struct get_adj_list_pmap<edge_owner_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename Graph::local_edge_descriptor
          local_edge_descriptor;
      public:
        typedef edge_owner_property_map<local_edge_descriptor> type;
        typedef type const_type;
      };
    };

    /**
     * The edge local property map maps from edge descriptors to
     * the local descriptor for that edge.
     */
    template<>
    struct get_adj_list_pmap<edge_local_t>
    {
      template<typename Graph>
      struct apply
      {
      private:
        typedef typename Graph::local_edge_descriptor
          local_edge_descriptor;
      public:
        typedef edge_local_property_map<local_edge_descriptor> type;
        typedef type const_type;
      };
    };
    //------------------------------------------------------------------------

    // Directed graphs do not have in edges, so this is a no-op
    template<typename Graph>
    inline void
    remove_in_edge(typename Graph::edge_descriptor, Graph&, directedS)
    { }

    //------------------------------------------------------------------------
    // Edge addition

    // Work around the fact that an adjacency_list with vecS vertex
    // storage automatically adds edges when the descriptor is
    // out-of-range.
    template <class Graph, class Config, class Base>
    inline std::pair<typename Config::edge_descriptor, bool>
    add_local_edge(typename Config::vertex_descriptor u,
                   typename Config::vertex_descriptor v,
                   const typename Config::edge_property_type& p,
                   vec_adj_list_impl<Graph, Config, Base>& g_)
    {
      adj_list_helper<Config, Base>& g = g_;
      return add_edge(u, v, p, g);
    }

    template <class Graph, class Config, class Base>
    inline std::pair<typename Config::edge_descriptor, bool>
    add_local_edge(typename Config::vertex_descriptor u,
                   typename Config::vertex_descriptor v,
                   const typename Config::edge_property_type& p,
                   boost::adj_list_impl<Graph, Config, Base>& g)
    {
      return add_edge(u, v, p, g);
    }

  } } // end namespace detail::parallel

  /**
   * Adjacency list traits for a distributed adjacency list. Contains
   * the vertex and edge descriptors, the directed-ness, and the
   * parallel edges typedefs.
   */
  template<typename OutEdgeListS, typename InVertexListS, 
	   typename InDistribution, typename DirectedS>
  struct adjacency_list_traits<OutEdgeListS,
                               distributedS<InVertexListS,
                                            InDistribution>,
                               DirectedS>
  {
  private:
    typedef typename mpl::if_<is_same<InVertexListS, defaultS>,
                              vecS,
                              InVertexListS>::type VertexListS;

    typedef adjacency_list_traits<OutEdgeListS, VertexListS, directedS>
      base_type;

  public:
    typedef typename base_type::vertex_descriptor local_vertex_descriptor;
    typedef typename base_type::edge_descriptor   local_edge_descriptor;

    typedef typename boost::mpl::if_<typename DirectedS::is_bidir_t,
      bidirectional_tag,
      typename boost::mpl::if_<typename DirectedS::is_directed_t,
        directed_tag, undirected_tag
      >::type
    >::type directed_category;

    typedef typename parallel_edge_traits<OutEdgeListS>::type
      edge_parallel_category;

    typedef detail::parallel::global_descriptor<local_vertex_descriptor>
      vertex_descriptor;

    typedef detail::parallel::edge_descriptor<local_edge_descriptor>
      edge_descriptor;
  };

#define PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS                                    \
  typename OutEdgeListS, typename InVertexListS, typename InDistribution,      \
  typename DirectedS, typename VertexProperty, typename EdgeProperty,          \
  typename GraphProperty, typename EdgeListS

#define PBGL_DISTRIB_ADJLIST_TYPE                                              \
  adjacency_list<OutEdgeListS,                                                 \
                 distributedS<InVertexListS, InDistribution>,                  \
                 DirectedS, VertexProperty, EdgeProperty, GraphProperty,       \
                 EdgeListS>

#define PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS_CONFIG                             \
  typename OutEdgeListS, typename InVertexListS, typename InDistribution,      \
  typename VertexProperty, typename EdgeProperty,  typename GraphProperty,     \
  typename EdgeListS
  
#define PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directed)                             \
  adjacency_list<OutEdgeListS,                                                 \
                 distributedS<InVertexListS, InDistribution>,                  \
                 directed, VertexProperty, EdgeProperty, GraphProperty,        \
                 EdgeListS>
                 
  /** A distributed adjacency list.
   *
   * This class template partial specialization defines a distributed
   * (or "partitioned") adjacency list graph. The distributed
   * adjacency list is similar to the standard Boost Graph Library
   * adjacency list, which stores a list of vertices and for each
   * verted the list of edges outgoing from the vertex (and, in some
   * cases, also the edges incoming to the vertex). The distributed
   * adjacency list differs in that it partitions the graph into
   * several subgraphs that are then divided among different
   * processors (or nodes within a cluster). The distributed adjacency
   * list attempts to maintain a high degree of compatibility with the
   * standard, non-distributed adjacency list.
   *
   * The graph is partitioned by vertex, with each processor storing
   * all of the required information for a particular subset of the
   * vertices, including vertex properties, outgoing edges, and (for
   * bidirectional graphs) incoming edges. This information is
   * accessible only on the processor that owns the vertex: for
   * instance, if processor 0 owns vertex @c v, no other processor can
   * directly access the properties of @c v or enumerate its outgoing
   * edges.
   *
   * Edges in a graph may be entirely local (connecting two local
   * vertices), but more often it is the case that edges are
   * non-local, meaning that the two vertices they connect reside in
   * different processes. Edge properties are stored with the
   * originating vertex for directed and bidirectional graphs, and are
   * therefore only accessible from the processor that owns the
   * originating vertex. Other processors may query the source and
   * target of the edge, but cannot access its properties. This is
   * particularly interesting when accessing the incoming edges of a
   * bidirectional graph, which are not guaranteed to be stored on the
   * processor that is able to perform the iteration. For undirected
   * graphs the situation is more complicated, since no vertex clearly
   * owns the edges: the list of edges incident to a vertex may
   * contain a mix of local and non-local edges.
   *
   * The distributed adjacency list is able to model several of the
   * existing Graph concepts. It models the Graph concept because it
   * exposes vertex and edge descriptors in the normal way; these
   * descriptors model the GlobalDescriptor concept (because they have
   * an owner and a local descriptor), and as such the distributed
   * adjacency list models the DistributedGraph concept. The adjacency
   * list also models the IncidenceGraph and AdjacencyGraph concepts,
   * although this is only true so long as the domain of the valid
   * expression arguments are restricted to vertices and edges stored
   * locally. Likewise, bidirectional and undirected distributed
   * adjacency lists model the BidirectionalGraph concept (vertex and
   * edge domains must be respectived) and the distributed adjacency
   * list models the MutableGraph concept (vertices and edges can only
   * be added or removed locally). T he distributed adjacency list
   * does not, however, model the VertexListGraph or EdgeListGraph
   * concepts, because we can not efficiently enumerate all vertices
   * or edges in the graph. Instead, the local subsets of vertices and
   * edges can be enumerated (with the same syntax): the distributed
   * adjacency list therefore models the DistributedVertexListGraph
   * and DistributedEdgeListGraph concepts, because concurrent
   * iteration over all of the vertices or edges stored on each
   * processor will visit each vertex or edge.
   *
   * The distributed adjacency list is distinguished from the
   * non-distributed version by the vertex list descriptor, which will
   * be @c distributedS<ProcessGroup,VertexListS>. Here,
   * the VertexListS type plays the same role as the VertexListS type
   * in the non-distributed adjacency list: it allows one to select
   * the data structure that will be used to store the local
   * vertices. The ProcessGroup type, on the other hand, is unique to
   * distributed data structures: it is the type that abstracts a
   * group of cooperating processes, and it used for process
   * identification, communication, and synchronization, among other
   * things. Different process group types represent different
   * communication mediums (e.g., MPI, PVM, TCP) or different models
   * of communication (LogP, CGM, BSP, synchronous, etc.). This
   * distributed adjacency list assumes a model based on non-blocking
   * sends.
   *
   * Distribution of vertices across different processors is
   * accomplished in two different ways. When initially constructing
   * the graph, the user may provide a distribution object (that
   * models the Distribution concept), which will determine the
   * distribution of vertices to each process. Additionally, the @c
   * add_vertex and @c add_edge operations add vertices or edges
   * stored on the local processor. For @c add_edge, this is
   * accomplished by requiring that the source vertex of the new edge
   * be local to the process executing @c add_edge.
   *
   * Internal properties of a distributed adjacency list are
   * accessible in the same manner as internal properties for a
   * non-distributed adjacency list for local vertices or
   * edges. Access to properties for remote vertices or edges occurs
   * with the same syntax, but involve communication with the owner of
   * the information: for more information, refer to class template
   * @ref distributed_property_map, which manages distributed
   * property maps. Note that the distributed property maps created
   * for internal properties determine their reduction operation via
   * the metafunction @ref property_reduce, which for the vast
   * majority of uses is correct behavior.
   *
   * Communication among the processes coordinating on a particular
   * distributed graph relies on non-blocking message passing along
   * with synchronization. Local portions of the distributed graph may
   * be modified concurrently, including the introduction of non-local
   * edges, but prior to accessing the graph it is recommended that
   * the @c synchronize free function be invoked on the graph to clear
   * up any pending interprocess communication and modifications. All
   * processes will then be released from the synchronization barrier
   * concurrently.
   *
   * \todo Determine precisely what we should do with nonlocal edges
   * in undirected graphs. Our parallelization of certain algorithms
   * relies on the ability to access edge property maps immediately
   * (e.g., edge_weight_t), so it may be necessary to duplicate the
   * edge properties in both processes (but then we need some form of
   * coherence protocol).
   *
   * \todo What does the user do if @c property_reduce doesn't do the
   * right thing?
   */
  template<typename OutEdgeListS, typename InVertexListS, 
	   typename InDistribution, typename DirectedS,
           typename VertexProperty, typename EdgeProperty, 
           typename GraphProperty, typename EdgeListS>
  class adjacency_list<OutEdgeListS,
                       distributedS<InVertexListS, 
                                    InDistribution>,
                       DirectedS, VertexProperty,
                       EdgeProperty, GraphProperty, EdgeListS>
  {
    typedef detail::parallel::adjacency_list_config<OutEdgeListS, InVertexListS, 
						    InDistribution, DirectedS, 
						    VertexProperty, EdgeProperty, 
						    GraphProperty, EdgeListS>
      config_type;
      
    typedef adjacency_list_traits<OutEdgeListS,
                                  distributedS<InVertexListS, 
                                               InDistribution>,
                                  DirectedS> 
      traits_type;

    typedef typename DirectedS::is_directed_t is_directed;

    typedef EdgeListS edge_list_selector;

  public:
    /// The type of the underlying adjacency list implementation
    typedef typename config_type::inherited inherited;
    typedef inherited base_type;

    /// The type of properties stored in the local subgraph
    /// Bidirectional graphs have an extra vertex property to store
    /// the incoming edges.
    typedef typename inherited::vertex_property_type
      base_vertex_property_type;

    /// The type of the distributed adjacency list (this type)
    typedef typename config_type::graph_type graph_type;

    /// Expose graph components and graph category
    typedef typename traits_type::local_vertex_descriptor
      local_vertex_descriptor;
    typedef typename traits_type::local_edge_descriptor
      local_edge_descriptor;
    typedef typename traits_type::vertex_descriptor vertex_descriptor;
    typedef typename traits_type::edge_descriptor edge_descriptor;

    typedef typename traits_type::directed_category directed_category;
    typedef typename inherited::edge_parallel_category
      edge_parallel_category;
    typedef typename inherited::graph_tag graph_tag;

    // Current implementation requires the ability to have parallel
    // edges in the underlying adjacency_list. Which processor each
    // edge refers to is attached as an internal property. TBD:
    // remove this restriction, which may require some rewriting.
    BOOST_STATIC_ASSERT((is_same<edge_parallel_category,
                                 allow_parallel_edge_tag>::value));

    /** Determine the graph traversal category.
     *
     * A directed distributed adjacency list models the Distributed
     * Graph, Incidence Graph, and Adjacency Graph
     * concepts. Bidirectional and undirected graphs also model the
     * Bidirectional Graph concept. Note that when modeling these
     * concepts the domains of certain operations (e.g., in_edges)
     * are restricted; see the distributed adjacency_list
     * documentation.
     */
    typedef typename boost::mpl::if_<
              is_same<DirectedS, directedS>,
              directed_distributed_adj_list_tag,
              typename boost::mpl::if_<is_same<DirectedS, bidirectionalS>,
                                       bidirectional_distributed_adj_list_tag,
                                       undirected_distributed_adj_list_tag>::type>
      ::type traversal_category;

    typedef typename inherited::degree_size_type degree_size_type;
    typedef typename inherited::vertices_size_type vertices_size_type;
    typedef typename inherited::edges_size_type edges_size_type;
    typedef VertexProperty vertex_property_type;
    typedef EdgeProperty edge_property_type;
    typedef typename inherited::graph_property_type graph_property_type;
    typedef typename inherited::vertex_bundled vertex_bundled;
    typedef typename inherited::edge_bundled edge_bundled;

    typedef typename container_gen<edge_list_selector, edge_descriptor>::type
      local_edge_list_type;

  private:
    typedef typename inherited::out_edge_iterator base_in_edge_iterator;

    typedef typename inherited::out_edge_iterator base_out_edge_iterator;
    typedef typename graph_traits<inherited>::edge_iterator
      base_edge_iterator;
    typedef typename inherited::edge_property_type base_edge_property_type;

    typedef typename local_edge_list_type::const_iterator
      undirected_edge_iterator;

    typedef InDistribution in_distribution_type;

  public:
    /// Iterator over the (local) vertices of the graph
    typedef transform_iterator<typename vertex_descriptor::generator,
                               typename inherited::vertex_iterator>
      vertex_iterator;

    /// Helper for out_edge_iterator
    typedef typename edge_descriptor::template out_generator<adjacency_list>
      out_edge_generator;

    /// Iterator over the outgoing edges of a vertex
    typedef transform_iterator<out_edge_generator,
                               typename inherited::out_edge_iterator>
      out_edge_iterator;

    /// Helper for in_edge_iterator
    typedef typename edge_descriptor::template in_generator<adjacency_list>
      in_edge_generator;

    /// Iterator over the incoming edges of a vertex
    typedef transform_iterator<in_edge_generator, base_in_edge_iterator>
      in_edge_iterator;

    /// Iterator over the neighbors of a vertex
    typedef boost::adjacency_iterator<
              adjacency_list, vertex_descriptor, out_edge_iterator,
              typename detail::iterator_traits<base_out_edge_iterator>
                         ::difference_type>
      adjacency_iterator;

    /// Iterator over the (local) edges in a graph
    typedef transform_iterator<out_edge_generator, base_edge_iterator>
      edge_iterator;

  public:
    /// How to refer to a process
    typedef typename amplusplus::transport::rank_type rank_type;

    /// Whether this graph is directed, undirected, or bidirectional
    typedef DirectedS directed_selector;

    /// the distribution type stored in the graph by default should be
    /// the variant_distribution type
    typedef typename mpl::if_<is_same<InDistribution, defaultS>,
			      boost::parallel::variant_distribution<vertices_size_type>,
			      InDistribution>::type base_distribution_type;

      typedef graph::distributed::shuffled_distribution<
          base_distribution_type> distribution_type;

    /// default_distribution_type is used as the distribution functor for the
    /// adjacency_list, it should be parallel::block by default
    typedef typename mpl::if_<is_same<InDistribution, defaultS>,
			      boost::parallel::block<local_vertex_descriptor>, base_distribution_type>::type
      default_distribution_type;

  private:
    // FIXME: the original adjacency_list contained this comment:
    //    Default copy constructor and copy assignment operators OK??? TBD
    // but the adj_list_impl contained these declarations:
    adjacency_list(const adjacency_list& other);
    adjacency_list& operator=(const adjacency_list& other);

  public:
    adjacency_list(amplusplus::transport& trans, 
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, default_distribution_type(trans, 0)),
        m_local_graph(GraphProperty()), 
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();
    }

    adjacency_list(amplusplus::transport& trans, 
                   const base_distribution_type& distribution,
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, distribution),
        m_local_graph(GraphProperty()), 
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();
    }

    adjacency_list(const GraphProperty& g,
                   amplusplus::transport& trans,
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, default_distribution_type(trans, 0)),
        m_local_graph(g), 
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();
    }

    adjacency_list(vertices_size_type n,
                   const GraphProperty& p,
                   amplusplus::transport& trans,
                   const base_distribution_type& distribution,
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, distribution),
        m_local_graph(distribution.block_size(trans.rank(), n), p),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
                                      get(vertex_index, base()));
    }

    adjacency_list(vertices_size_type n,
                   amplusplus::transport& trans,
                   const base_distribution_type& distribution,
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, distribution),
        m_local_graph(distribution.block_size(trans.rank(), n), GraphProperty()),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
                                      get(vertex_index, base()));
    }

    adjacency_list(vertices_size_type n,
                   const GraphProperty& p,
                   amplusplus::transport& trans,
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, default_distribution_type(trans, n)),
        m_local_graph(default_distribution_type(trans, n).block_size(trans.rank(), n), p),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
                                      get(vertex_index, base()));
    }

    adjacency_list(vertices_size_type n,
                   amplusplus::transport& trans,
		   int coalescing_size = 1 << 12)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, default_distribution_type(trans, n)),
        m_local_graph(default_distribution_type(trans, n).block_size(trans.rank(), n), 
                      GraphProperty()),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
                                      get(vertex_index, base()));
    }

    /*
     * We assume that every processor sees the same list of edges, so
     * they skip over any that don't originate from themselves. This
     * means that programs switching between a local and a distributed
     * graph will keep the same semantics.
     */
    template <class EdgeIterator>
    adjacency_list(EdgeIterator first, EdgeIterator last,
                   vertices_size_type n,
                   amplusplus::transport& trans, 
		   int coalescing_size = 1 << 12,
                   const GraphProperty& p = GraphProperty())
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, default_distribution_type(trans, n)),
        m_local_graph(default_distribution_type(trans, n).block_size(trans.rank(), n), p),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      {
	amplusplus::scoped_epoch epoch(transport_);
	typedef typename config_type::VertexListS vertex_list_selector;
	initialize(first, last, n, this->distribution(), vertex_list_selector());
	detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
                                      get(vertex_index, base()));
      }
    }

    template <class EdgeIterator, class EdgePropertyIterator>
    adjacency_list(EdgeIterator first, EdgeIterator last,
                   EdgePropertyIterator ep_iter,
                   vertices_size_type n,
                   amplusplus::transport& trans,
		   int coalescing_size = 1 << 12,
                   const GraphProperty& p = GraphProperty())
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, default_distribution_type(trans, n)),
        m_local_graph(default_distribution_type(trans, n).block_size(trans.rank(), n), p),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      {
	amplusplus::scoped_epoch epoch(transport_);
	typedef typename config_type::VertexListS vertex_list_selector;
	initialize(first, last, ep_iter, n, this->distribution(),
		   vertex_list_selector());
	detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
                                      get(vertex_index, base()));
      }
    }

    template <class EdgeIterator>
    adjacency_list(EdgeIterator first, EdgeIterator last,
                   vertices_size_type n,
                   amplusplus::transport& trans,
                   const base_distribution_type& distribution,
		   int coalescing_size = 1 << 12,
                   const GraphProperty& p = GraphProperty())
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, distribution),
        m_local_graph(distribution.block_size(trans.rank(), n), p),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      {
	amplusplus::scoped_epoch epoch(transport_);
	typedef typename config_type::VertexListS vertex_list_selector;
	initialize(first, last, n, this->distribution(), vertex_list_selector());
	detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
							  get(vertex_index, base()));
      }
    }

    template <class EdgeIterator, class EdgePropertyIterator>
    adjacency_list(EdgeIterator first, EdgeIterator last,
                   EdgePropertyIterator ep_iter,
                   vertices_size_type n,
                   amplusplus::transport& trans,
                   const base_distribution_type& distribution,
		   int coalescing_size = 1 << 12,
                   const GraphProperty& p = GraphProperty())
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_add_edge_data>(),
	   amplusplus::register_mpi_datatype<msg_add_edge_with_property_data>(), 
	   amplusplus::register_mpi_datatype<msg_remove_edge_data>(), 0)),
        distribution_(trans, distribution),
        m_local_graph(distribution.block_size(trans.rank(), n), p),
        transport_(trans),
	add_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	add_edge_with_property_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans),
	remove_edge_msg(amplusplus::basic_coalesced_message_type_gen(coalescing_size), trans)
    {
      setup_messages();

      {
	amplusplus::scoped_epoch epoch(transport_);
	typedef typename config_type::VertexListS vertex_list_selector;
	initialize(first, last, ep_iter, n, distribution,
		   vertex_list_selector());
	detail::parallel::maybe_initialize_vertex_indices(vertices(base()),
							  get(vertex_index, base()));
      }
    }

    void clear()
    {
      base().clear();
    }

    static vertex_descriptor null_vertex()
    {
      return vertex_descriptor(rank_type(0), inherited::null_vertex());
    }

    inherited&       base()       { return m_local_graph; }
    const inherited& base() const { return m_local_graph; }

    rank_type processor() const { return transport_.rank(); }

    amplusplus::transport&       transport() { return transport_; }
    const amplusplus::transport& transport() const { return transport_; }

    distribution_type&       distribution()       { return this->distribution_; }
    const distribution_type& distribution() const { return this->distribution_; }

    // Redistribute the vertices of the graph by placing each vertex
    // v on the processor get(vertex_to_processor, v).
    template<typename VertexProcessorMap>
    void redistribute(VertexProcessorMap vertex_to_processor) { assert(false); }

    // Directly access a vertex or edge bundle
    vertex_bundled& operator[](vertex_descriptor v)
    {
      assert(v.owner == processor());
      return base()[v.local];
    }
    
    const vertex_bundled& operator[](vertex_descriptor v) const
    {
      assert(v.owner == processor());
      return base()[v.local];
    }
    
    edge_bundled& operator[](edge_descriptor e)
    {
      assert(e.owner() == processor());
      return base()[e.local];
    }
    
    const edge_bundled& operator[](edge_descriptor e) const
    {
      assert(e.owner() == processor());
      return base()[e.local];
    }

    template<typename OStreamConstructibleArchive>
    void save(std::string const& filename) const;

    template<typename IStreamConstructibleArchive>
    void load(std::string const& filename);

  private:
    // Request vertex->processor mapping for neighbors <does nothing>
    template<typename VertexProcessorMap>
    void
    request_in_neighbors(vertex_descriptor,
                         VertexProcessorMap,
                         directedS) { }

    // Clear the list of in-edges, but don't tell the remote processor
    void clear_in_edges_local(vertex_descriptor v, directedS) {}

    // Remove in-edges that have migrated <does nothing>
    template<typename VertexProcessorMap>
    void
    remove_migrated_in_edges(vertex_descriptor,
                             VertexProcessorMap,
                             directedS) { }

    // Initialize the graph with the given edge list and vertex
    // distribution. This variation works only when
    // VertexListS=vecS, and we know how to create remote vertex
    // descriptors based solely on the distribution.
    template<typename EdgeIterator>
    void
    initialize(EdgeIterator first, EdgeIterator last,
               vertices_size_type, const base_distribution_type& distribution, 
               vecS);

    // Initialize the graph with the given edge list, edge
    // properties, and vertex distribution. This variation works
    // only when VertexListS=vecS, and we know how to create remote
    // vertex descriptors based solely on the distribution.
    template<typename EdgeIterator, typename EdgePropertyIterator>
    void
    initialize(EdgeIterator first, EdgeIterator last,
               EdgePropertyIterator ep_iter,
               vertices_size_type, const base_distribution_type& distribution, 
               vecS);

    // Initialize the graph with the given edge list, edge
    // properties, and vertex distribution.
    template<typename EdgeIterator, typename EdgePropertyIterator,
             typename VertexListS>
    void
    initialize(EdgeIterator first, EdgeIterator last,
               EdgePropertyIterator ep_iter,
               vertices_size_type n, 
               const base_distribution_type& distribution,
               VertexListS);

    // Initialize the graph with the given edge list and vertex
    // distribution. This is nearly identical to the one below it,
    // for which I should be flogged. However, this version does use
    // slightly less memory than the version that accepts an edge
    // property iterator.
    template<typename EdgeIterator, typename VertexListS>
    void
    initialize(EdgeIterator first, EdgeIterator last,
               vertices_size_type n, 
               const base_distribution_type& distribution,
               VertexListS);

  public:
    //---------------------------------------------------------------------
    // Build a vertex property instance for the underlying adjacency
    // list from the given property instance of the type exposed to
    // the user.
    base_vertex_property_type 
    build_vertex_property(const vertex_property_type& p)
    { return build_vertex_property(p, directed_selector()); }

    base_vertex_property_type
    build_vertex_property(const vertex_property_type& p, directedS)
    {
      return base_vertex_property_type(p);
    }

    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Build an edge property instance for the underlying adjacency
    // list from the given property instance of the type exposed to
    // the user.
    base_edge_property_type build_edge_property(const edge_property_type& p)
    { return build_edge_property(p, directed_selector()); }

    base_edge_property_type
    build_edge_property(const edge_property_type& p, directedS)
    {
      return base_edge_property_type(0, p);
    }

  private:
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Opposite of above.
    edge_property_type split_edge_property(const base_edge_property_type& p)
    { return split_edge_property(p, directed_selector()); }

    edge_property_type
    split_edge_property(const base_edge_property_type& p, directedS)
    {
      return p.m_base;
    }

    edge_property_type
    split_edge_property(const base_edge_property_type& p, bidirectionalS)
    {
      return p.m_base;
    }

    edge_property_type
    split_edge_property(const base_edge_property_type& p, undirectedS)
    {
      return p.m_base.m_base;
    }
    //---------------------------------------------------------------------

    /** The set of messages that can be transmitted and received by
     *  a distributed adjacency list. This list will eventually be
     *  exhaustive, but is currently quite limited.
     */

    struct add_edge_handler;
    friend struct add_edge_handler;

    struct add_edge_with_property_handler;
    friend struct add_edge_with_property_handler;

    struct remove_edge_handler;
    friend struct remove_edge_handler;

    /// Message data structures

  public: /* Needed so free functions add_vertex(), add_edge(), etc. can send messages */

    typedef std::pair<local_vertex_descriptor, vertex_descriptor> msg_add_edge_data;
    typedef boost::tuple<local_vertex_descriptor, vertex_descriptor, edge_property_type>
      msg_add_edge_with_property_data;
    typedef edge_descriptor msg_remove_edge_data;

    // NGE: If vertex_property_type is 'no_property' then MPI datatype construction fails
    //      Adding vertices remotely is UNUSED at the moment and the message type is never 
    //      registered thus avoiding this problem for the moment.  The fix would be to send 
    //      a different message type when the vertex property is no_property

  private:
    /// Message types
    typedef amplusplus::basic_coalesced_message_type<msg_add_edge_data,
						     add_edge_handler> 
      add_edge_message_type;

    typedef amplusplus::basic_coalesced_message_type<msg_add_edge_with_property_data,
						     add_edge_with_property_handler> 
      add_edge_with_property_message_type;

    typedef amplusplus::basic_coalesced_message_type<msg_remove_edge_data,
						     remove_edge_handler> 
      remove_edge_message_type;

    /// Process incoming messages.
    void setup_messages();

  public:
    void send_remove_edge_request(edge_descriptor e)
    {
      rank_type dest = e.target_processor;
      if (e.target_processor == transport_.rank())
        dest = e.source_processor;
      remove_edge_msg.send(e, dest);
    }

  private:
    const int dummy_first_member_for_init_order; // Unused

    /// The distribution we will use to map names to processors
    distribution_type distribution_;

    /// The local subgraph
    inherited m_local_graph;

    /// The process group through which this distributed graph
    /// communicates.
    amplusplus::transport& transport_;

  public: /* Needed so free functions add_vertex(), add_edge(), etc. can send messages) */
    /// Messages
    add_edge_message_type                 add_edge_msg;
    add_edge_with_property_message_type   add_edge_with_property_msg;
    remove_edge_message_type              remove_edge_msg;
  };

  /**
   * Returns the set of vertices local to this processor. Note that
   * although this routine matches a valid expression of a
   * VertexListGraph, it does not meet the semantic requirements of
   * VertexListGraph because it returns only local vertices (not all
   * vertices).
   */
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  std::pair<typename PBGL_DISTRIB_ADJLIST_TYPE
                       ::vertex_iterator,
            typename PBGL_DISTRIB_ADJLIST_TYPE
                       ::vertex_iterator>
  vertices(const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE
      ::vertex_descriptor Vertex;

    typedef typename Vertex::generator generator;

    return std::make_pair(make_transform_iterator(vertices(g.base()).first,
                                                  generator(g.processor())),
                          make_transform_iterator(vertices(g.base()).second,
                                                  generator(g.processor())));
  }

  /**
   * Returns the number of vertices local to this processor. Note that
   * although this routine matches a valid expression of a
   * VertexListGraph, it does not meet the semantic requirements of
   * VertexListGraph because it returns only a count of local vertices
   * (not all vertices).
   */
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename PBGL_DISTRIB_ADJLIST_TYPE
             ::vertices_size_type
  num_vertices(const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    return num_vertices(g.base());
  }

  /***************************************************************************
   * Implementation of Incidence Graph concept
   ***************************************************************************/
  /**
   * Returns the source of edge @param e in @param g.
   */
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename Edge>
  typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor
  source(const detail::parallel::edge_descriptor<Edge>& e,
         const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE
      ::vertex_descriptor Vertex;
    return Vertex(e.source_processor, source(e.local, g.base()));
  }

  /**
   * Returns the target of edge @param e in @param g.
   */
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename Edge>
  typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor
  target(const detail::parallel::edge_descriptor<Edge>& e,
         const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE
      ::vertex_descriptor Vertex;
    return Vertex(e.target_processor, target(e.local, g.base()));
  }

  /**
   * Return the set of edges outgoing from a particular vertex. The
   * vertex @param v must be local to the processor executing this
   * routine.
   */
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  std::pair<typename PBGL_DISTRIB_ADJLIST_TYPE::out_edge_iterator,
            typename PBGL_DISTRIB_ADJLIST_TYPE::out_edge_iterator>
  out_edges(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
            const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    assert(v.owner == g.processor());

    typedef PBGL_DISTRIB_ADJLIST_TYPE impl;
    typedef typename impl::out_edge_generator generator;

    return std::make_pair(
             make_transform_iterator(out_edges(v.local, g.base()).first,
                                     generator(g)),
             make_transform_iterator(out_edges(v.local, g.base()).second,
                                     generator(g)));
  }

  /**
   * Return the number of edges outgoing from a particular vertex. The
   * vertex @param v must be local to the processor executing this
   * routine.
   */
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename PBGL_DISTRIB_ADJLIST_TYPE::degree_size_type
  out_degree(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
             const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    assert(v.owner == g.processor());

    return out_degree(v.local, g.base());
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename PBGL_DISTRIB_ADJLIST_TYPE::edges_size_type
  num_edges(const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    return num_edges(g.base());
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  std::pair<
    typename PBGL_DISTRIB_ADJLIST_TYPE::edge_iterator,
    typename PBGL_DISTRIB_ADJLIST_TYPE::edge_iterator>
  edges(const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE impl;
    typedef typename impl::out_edge_generator generator;

    return std::make_pair(make_transform_iterator(edges(g.base()).first,
                                                  generator(g)),
                          make_transform_iterator(edges(g.base()).second,
                                                  generator(g)));
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  inline
  typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor
  vertex(typename PBGL_DISTRIB_ADJLIST_TYPE::vertices_size_type n,
         const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor 
      vertex_descriptor;

    return vertex_descriptor(g.distribution()(n), g.distribution().local(n));
  }

  /***************************************************************************
   * Access to particular edges
   ***************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS_CONFIG>
  std::pair<
    typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)::edge_descriptor,
    bool
  >
  edge(typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)::vertex_descriptor u,
       typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)::vertex_descriptor v,
       const PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)
                       ::edge_descriptor edge_descriptor;

    // For directed graphs, u must be local
    assert(u.owner == g.transport().rank());

    typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)
        ::out_edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
      if (target(*ei, g) == v) return std::make_pair(*ei, true);
    }
    return std::make_pair(edge_descriptor(), false);
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  std::pair<
    typename PBGL_DISTRIB_ADJLIST_TYPE::edge_descriptor,
    bool
  >
  edge(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
       typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
       const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE
                       ::edge_descriptor edge_descriptor;

    // For bidirectional and undirected graphs, u must be local or v
    // must be local
    if (u.owner == g.transport().rank()) {
      typename PBGL_DISTRIB_ADJLIST_TYPE::out_edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
        if (target(*ei, g) == v) return std::make_pair(*ei, true);
      }
      return std::make_pair(edge_descriptor(), false);
    } else if (v.owner == g.transport().rank()) {
      typename PBGL_DISTRIB_ADJLIST_TYPE::in_edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = in_edges(v, g); ei != ei_end; ++ei) {
        if (source(*ei, g) == u) return std::make_pair(*ei, true);
      }
      return std::make_pair(edge_descriptor(), false);
    } else {
      assert(false);
      exit(1);
    }
  }

#if 0
  // TBD: not yet supported
  std::pair<out_edge_iterator, out_edge_iterator>
  edge_range(vertex_descriptor u, vertex_descriptor v,
             const adjacency_list& g);
#endif

  /***************************************************************************
   * Implementation of Adjacency Graph concept
   ***************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  std::pair<typename PBGL_DISTRIB_ADJLIST_TYPE::adjacency_iterator,
            typename PBGL_DISTRIB_ADJLIST_TYPE::adjacency_iterator>
  adjacent_vertices(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
                    const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE::adjacency_iterator iter;
    return std::make_pair(iter(out_edges(v, g).first, &g),
                          iter(out_edges(v, g).second, &g));
  }

  /***************************************************************************
   * Implementation of Mutable Graph concept
   ***************************************************************************/

  /************************************************************************
   * add_edge
   ************************************************************************/

  // NGE: Having add edge return a descriptor requires communicating
  //      remote edge insertions back to the requesting processor,
  //      which is unnecessary in most cases.  The bookkeeping to
  //      allow the commit-on-dereference implementation in the
  //      original PBGL is a source of significant overhead even if it
  //      is never used.

  //      Undirected graphs would require similar techniques, along
  //      with a method for matching edges and ensuring property
  //      consistency.  However it's my opinion that this is better
  //      implemented as a separate graph class due to the significant
  //      additional complexity.

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  add_edge(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
           typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
           PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typename PBGL_DISTRIB_ADJLIST_TYPE::edge_property_type e;
    add_edge(u, v, e, g);
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  add_edge(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
           typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
           typename PBGL_DISTRIB_ADJLIST_TYPE::edge_property_type const& p,
           PBGL_DISTRIB_ADJLIST_TYPE& g)
  {

    if (u.owner == g.processor()) {

      // Add the edge to the local part of the graph
      std::pair<typename PBGL_DISTRIB_ADJLIST_TYPE::local_edge_descriptor, bool> 
	inserted = detail::parallel::add_local_edge(u.local, v.local,
						    g.build_edge_property(p), 
						    g.base());

      if (inserted.second)
	// Keep track of the owner of the target
	put(edge_target_processor_id, g.base(), inserted.first, 
	    v.owner);
    } else
      // Request that the remote processor add an edge and, but
      // don't wait for a reply.
      g.add_edge_msg.send(std::make_pair(u.local, v), u.owner);
  }

  /************************************************************************
   *
   * remove_edge
   *
   ************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  remove_edge(typename PBGL_DISTRIB_ADJLIST_TYPE::edge_descriptor e,
              PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    assert(source(e, g).owner == g.processor()
           || target(e, g).owner == g.processor());

    if (target(e, g).owner == g.processor())
      detail::parallel::remove_in_edge(e, g, DirectedS());
    if (source(e, g).owner == g.processor())
      remove_edge(e.local, g.base());

    if (source(e, g).owner != g.processor()
        || (target(e, g).owner != g.processor()
            && !(is_same<DirectedS, directedS>::value))) {
      g.send_remove_edge_request(e);
    }
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  remove_edge(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
              typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v,
              PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE
                       ::vertex_descriptor vertex_descriptor;
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE
                       ::edge_descriptor edge_descriptor;
    std::pair<edge_descriptor, bool> the_edge = edge(u, v, g);
    if (the_edge.second) remove_edge(the_edge.first, g);
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  inline void
  remove_edge(typename PBGL_DISTRIB_ADJLIST_TYPE::out_edge_iterator ei,
              PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    remove_edge(*ei, g);
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS_CONFIG>
  inline void
  remove_edge(typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)
                ::out_edge_iterator ei,
              PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)& g)
  {
    assert(source(*ei, g).owner == g.processor());
    remove_edge(ei->local, g.base());
  }

  /************************************************************************
   *
   * remove_out_edge_if
   *
   ************************************************************************/
  namespace parallel { namespace detail {
    /**
     * Function object that applies the underlying predicate to
     * determine if an out-edge should be removed. If so, either
     * removes the incoming edge (if it is stored locally) or sends a
     * message to the owner of the target requesting that it remove
     * the edge.
     */
    template<typename Graph, typename Predicate>
    struct remove_out_edge_predicate
    {
      typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
      typedef typename Graph::local_edge_descriptor         argument_type;
      typedef typename Graph::directed_selector             directed_selector;
      typedef bool result_type;

      remove_out_edge_predicate(Graph& g, Predicate& predicate)
        : g(g), predicate(predicate) { }

      bool operator()(const argument_type& le)
      {
        typedef typename edge_descriptor::template out_generator<Graph>
          generator;

        edge_descriptor e = generator(g)(le);

        if (predicate(e)) {
          if (source(e, g).owner != target(e, g).owner
              && !(is_same<directed_selector, directedS>::value))
            g.send_remove_edge_request(e);
          else
            ::boost::detail::parallel::remove_in_edge(e, g,
                                                      directed_selector());
          return true;
        } else return false;
      }

    private:
      Graph& g;
      Predicate predicate;
    };
  } } // end namespace parallel::detail

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename Predicate>
  inline void
  remove_out_edge_if
     (typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
      Predicate predicate,
      PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE Graph;
    typedef parallel::detail::remove_out_edge_predicate<Graph, Predicate>
      Pred;

    assert(u.owner == g.processor());
    remove_out_edge_if(u.local, Pred(g, predicate), g.base());
  }

  /************************************************************************
   *
   * remove_edge_if
   *
   ************************************************************************/
  namespace parallel { namespace detail {
    /**
     * Function object that applies the underlying predicate to
     * determine if a directed edge can be removed. This only applies
     * to directed graphs.
     */
    template<typename Graph, typename Predicate>
    struct remove_directed_edge_predicate
    {
      typedef typename Graph::local_edge_descriptor argument_type;
      typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
      typedef bool result_type;

      remove_directed_edge_predicate(Graph& g, const Predicate& predicate)
        : g(g), predicate(predicate) { }

      bool operator()(const argument_type& le)
      {
        typedef typename edge_descriptor::template out_generator<Graph>
          generator;

        edge_descriptor e = generator(g)(le);
        return predicate(e);
      }

    private:
      Graph& g;
      Predicate predicate;
    };
  } } // end namespace parallel::detail

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS_CONFIG, typename Predicate>
  inline void
  remove_edge_if(Predicate predicate, 
                 PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS) Graph;
    typedef parallel::detail::remove_directed_edge_predicate<Graph,
                                                             Predicate> Pred;
    remove_edge_if(Pred(g, predicate), g.base());
  }

  /************************************************************************
   *
   * clear_vertex
   *
   ************************************************************************/
  namespace parallel { namespace detail {
    struct always_true
    {
      typedef bool result_type;

      template<typename T> bool operator()(const T&) const { return true; }
    };
  } } // end namespace parallel::detail

  /************************************************************************
   *
   * clear_out_edges
   *
   ************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS_CONFIG>
  void
  clear_out_edges
    (typename PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)::vertex_descriptor u,
      PBGL_DISTRIB_ADJLIST_TYPE_CONFIG(directedS)& g)
  {
    assert(u.owner == g.processor());
    clear_out_edges(u.local, g.base());
  }

  /************************************************************************
   *
   * add_vertex
   *
   ************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  add_vertex(PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_property_type p;
    add_vertex(p, g);
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  add_vertex(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_property_type const& p,
             PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    // Add vertex locally since we have no owner_by_property from named_graph
    add_vertex(g.build_vertex_property(p), g.base());
  }

  /************************************************************************
   *
   * remove_vertex
   *
   ************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  void
  remove_vertex(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
                PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef typename PBGL_DISTRIB_ADJLIST_TYPE::graph_type graph_type;
    assert(u.owner == g.processor());
    g.distribution().clear();
    remove_vertex(u.local, g.base());
  }

  /***************************************************************************
   * Implementation of Property Graph concept
   ***************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename Property>
  struct property_map<PBGL_DISTRIB_ADJLIST_TYPE, Property>
  : detail::parallel::get_adj_list_pmap<Property>
      ::template apply<PBGL_DISTRIB_ADJLIST_TYPE>
  { };

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename Property>
  struct property_map<PBGL_DISTRIB_ADJLIST_TYPE const, Property>
          : boost::detail::parallel::get_adj_list_pmap<Property>
// FIXME: in the original code the following was not const
      ::template apply<PBGL_DISTRIB_ADJLIST_TYPE const>
  { };

  template<typename Property, PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, Property>::type
  get(Property p, PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE Graph;
    typedef typename property_map<Graph, Property>::type result_type;
    typedef typename property_traits<result_type>::value_type value_type;
    typedef typename property_reduce<Property>::template apply<value_type>
      reduce;

    typedef typename property_traits<result_type>::key_type descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename mpl::if_<is_same<descriptor, vertex_descriptor>,
                              vertex_global_t, edge_global_t>::type
      global_map_t;

    return result_type(g.transport(), get(global_map_t(), g),
                       get(p, g.base()), reduce());
  }

  template<typename Property, PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, Property>::const_type
  get(Property p, const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE Graph;
    typedef typename property_map<Graph, Property>::const_type result_type;
    typedef typename property_traits<result_type>::value_type value_type;
    typedef typename property_reduce<Property>::template apply<value_type>
      reduce;

    typedef typename property_traits<result_type>::key_type descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename mpl::if_<is_same<descriptor, vertex_descriptor>,
                              vertex_global_t, edge_global_t>::type
      global_map_t;

    return result_type(g.transport(), get(global_map_t(), g),
                       get(p, g.base()), reduce());
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_local_index_t>::type
  get(vertex_local_index_t, PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    return get(vertex_local_index, g.base());
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE,
                        vertex_local_index_t>::const_type
  get(vertex_local_index_t, const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    return get(vertex_local_index, g.base());
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_global_t>::const_type
  get(vertex_global_t, const PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       vertex_global_t>::const_type result_type;
    return result_type();
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_global_t>::const_type
  get(vertex_global_t, PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       vertex_global_t>::const_type result_type;
    return result_type();
  }

  /// Retrieve a property map mapping from a vertex descriptor to its
  /// owner.
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_owner_t>::type
  get(vertex_owner_t, PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       vertex_owner_t>::type result_type;
    return result_type();
  }

  /// Retrieve a property map mapping from a vertex descriptor to its
  /// owner.
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_owner_t>::const_type
  get(vertex_owner_t, const PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       vertex_owner_t>::const_type result_type;
    return result_type();
  }

  /// Retrieve the owner of a vertex
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  inline rank_type
  get(vertex_owner_t, PBGL_DISTRIB_ADJLIST_TYPE&,
      typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v)
  {
    return v.owner;
  }

  /// Retrieve the owner of a vertex
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  inline rank_type
  get(vertex_owner_t, const PBGL_DISTRIB_ADJLIST_TYPE&,
      typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v)
  {
    return v.owner;
  }

  /// Retrieve a property map that maps from a vertex descriptor to
  /// its local descriptor.
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_local_t>::type
  get(vertex_local_t, PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       vertex_local_t>::type result_type;
    return result_type();
  }

  /// Retrieve a property map that maps from a vertex descriptor to
  /// its local descriptor.
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_local_t>::const_type
  get(vertex_local_t, const PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       vertex_local_t>::const_type result_type;
    return result_type();
  }

  /// Retrieve the local descriptor of a vertex
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  inline typename PBGL_DISTRIB_ADJLIST_TYPE::local_vertex_descriptor
  get(vertex_local_t, PBGL_DISTRIB_ADJLIST_TYPE&,
      typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v)
  {
    return v.local;
  }

  /// Retrieve the local descriptor of a vertex
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  inline typename PBGL_DISTRIB_ADJLIST_TYPE::local_vertex_descriptor
  get(vertex_local_t, const PBGL_DISTRIB_ADJLIST_TYPE&,
      typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor v)
  {
    return v.local;
  }


  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, edge_global_t>::const_type
  get(edge_global_t, const PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       edge_global_t>::const_type result_type;
    return result_type();
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, edge_global_t>::const_type
  get(edge_global_t, PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       edge_global_t>::const_type result_type;
    return result_type();
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, edge_owner_t>::type
  get(edge_owner_t, PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       edge_owner_t>::type result_type;
    return result_type();
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, edge_owner_t>::const_type
  get(edge_owner_t, const PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       edge_owner_t>::const_type result_type;
    return result_type();
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, edge_local_t>::type
  get(edge_local_t, PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       edge_local_t>::type result_type;
    return result_type();
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, edge_local_t>::const_type
  get(edge_local_t, const PBGL_DISTRIB_ADJLIST_TYPE& )
  {
    typedef typename property_map<
                       PBGL_DISTRIB_ADJLIST_TYPE,
                       edge_local_t>::const_type result_type;
    return result_type();
  }

  template<typename Property, PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS,
           typename Key>
  inline
  typename property_traits<typename property_map<
                PBGL_DISTRIB_ADJLIST_TYPE, Property>::const_type
           >::value_type
  get(Property p, const PBGL_DISTRIB_ADJLIST_TYPE& g, const Key& key)
  {
    if (owner(key) == g.transport().rank())
      return get(p, g.base(), local(key));
    else
      assert(false);
  }

  template<typename Property, PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS,
           typename Key, typename Value>
  void
  put(Property p, PBGL_DISTRIB_ADJLIST_TYPE& g, const Key& key, const Value& v)
  {
    if (owner(key) == g.transport().rank())
      put(p, g.base(), local(key), v);
    else
      assert(false);
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_index_t>::type
  get(vertex_index_t vi, PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE graph_type;
    typedef typename property_map<graph_type, vertex_index_t>::type
      result_type;
    return result_type(g.transport(), get(vertex_global, g),
                       get(vi, g.base()));
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, vertex_index_t>::const_type
  get(vertex_index_t vi, const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE graph_type;
    typedef typename property_map<graph_type, vertex_index_t>::const_type
      result_type;
    return result_type(g.transport(), get(vertex_global, g),
                       get(vi, g.base()));
  }

#if 0
  /***************************************************************************
   * Implementation of bundled properties
   ***************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename T, typename Bundle>
  struct property_map<PBGL_DISTRIB_ADJLIST_TYPE, T Bundle::*>
    : detail::parallel::get_adj_list_pmap<T Bundle::*>
      ::template apply<PBGL_DISTRIB_ADJLIST_TYPE>
  { };

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename T, typename Bundle>
  struct property_map<PBGL_DISTRIB_ADJLIST_TYPE const, T Bundle::*>
    : detail::parallel::get_adj_list_pmap<T Bundle::*>
      ::template apply<PBGL_DISTRIB_ADJLIST_TYPE const>
  { };

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename T, typename Bundle>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, T Bundle::*>::type
  get(T Bundle::* p, PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE Graph;
    typedef typename property_map<Graph, T Bundle::*>::type result_type;
    typedef typename property_traits<result_type>::value_type value_type;
    typedef typename property_reduce<T Bundle::*>::template apply<value_type>
      reduce;

    typedef typename property_traits<result_type>::key_type descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename mpl::if_<is_same<descriptor, vertex_descriptor>,
                              vertex_global_t, edge_global_t>::type
      global_map_t;

    return result_type(g.transport(), get(global_map_t(), g),
                       get(p, g.base()), reduce());
  }

  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS, typename T, typename Bundle>
  typename property_map<PBGL_DISTRIB_ADJLIST_TYPE, T Bundle::*>::const_type
  get(T Bundle::* p, const PBGL_DISTRIB_ADJLIST_TYPE& g)
  {
    typedef PBGL_DISTRIB_ADJLIST_TYPE Graph;
    typedef typename property_map<Graph, T Bundle::*>::const_type result_type;
    typedef typename property_traits<result_type>::value_type value_type;
    typedef typename property_reduce<T Bundle::*>::template apply<value_type>
      reduce;

    typedef typename property_traits<result_type>::key_type descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename mpl::if_<is_same<descriptor, vertex_descriptor>,
                              vertex_global_t, edge_global_t>::type
      global_map_t;

    return result_type(g.transport(), get(global_map_t(), g),
                       get(p, g.base()), reduce());
  }
#endif

  /***************************************************************************
   * Implementation of DistributedGraph concept
   ***************************************************************************/
  template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
  amplusplus::transport
  process_group(const PBGL_DISTRIB_ADJLIST_TYPE& g)
  { return g.transport(); }

  // Hash function for global descriptors
  template<typename LocalDescriptor>
  struct hash<detail::parallel::global_descriptor<LocalDescriptor> >
  {
    typedef detail::parallel::global_descriptor<LocalDescriptor> argument_type;
    std::size_t operator()(argument_type const& x) const
    {
      std::size_t hash = hash_value(x.owner);
      hash_combine(hash, x.local);
      return hash;
    }
  };

  // Hash function for parallel edge descriptors
  template<typename Edge>
  struct hash<detail::parallel::edge_descriptor<Edge> >
  {
    typedef detail::parallel::edge_descriptor<Edge> argument_type;

    std::size_t operator()(argument_type const& x) const
    {
      std::size_t hash = hash_value(x.owner());
      hash_combine(hash, x.local);
      return hash;
    }
  };

} // end namespace boost

#include <boost/graph/distributed/adjlist/handlers.hpp>
#include <boost/graph/distributed/adjlist/initialize.hpp>

namespace amplusplus {
  // Global descriptor
  template <typename A>
  struct make_mpi_datatype<boost::detail::parallel::global_descriptor<A> > 
    : make_mpi_datatype_base 
  {
    make_mpi_datatype<std::size_t> dt1;
    make_mpi_datatype<A> dt2;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1(), dt2() {
      int blocklengths[2] = {1, 1};
      MPI_Aint displacements[2];
      char dummy;
      boost::detail::parallel::global_descriptor<A> *test_object
	= (boost::detail::parallel::global_descriptor<A>*)(&dummy);
      MPI_Aint test_object_ptr;
      MPI_Get_address(test_object, &test_object_ptr);
      MPI_Get_address(&test_object->owner, &displacements[0]);
      MPI_Get_address(&test_object->local, &displacements[1]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      MPI_Datatype types[2] = {dt1.get(), dt2.get()};
      MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };

  template <typename A>
  struct dummy_value<boost::detail::parallel::global_descriptor<A>, void> {
    boost::detail::parallel::global_descriptor<A>
    operator()(size_t) const { return boost::detail::parallel::global_descriptor<A>(0, std::numeric_limits<A>::max()); }
  };

  // Distributed edge descriptors
  template <typename A>
  struct make_mpi_datatype<boost::detail::parallel::edge_descriptor<A> > 
    : make_mpi_datatype_base 
  {
    make_mpi_datatype<std::size_t> dt1;
    make_mpi_datatype<std::size_t> dt2;
    make_mpi_datatype<A> dt3;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1(), dt2(), dt3() {
      int blocklengths[3] = {1, 1, 1};
      MPI_Aint displacements[3];
      char dummy;
      boost::detail::parallel::edge_descriptor<A> *test_object 
	= (boost::detail::parallel::edge_descriptor<A>*)(&dummy);
      MPI_Aint test_object_ptr;
      MPI_Get_address(test_object, &test_object_ptr);
      MPI_Get_address(&test_object->source_processor, &displacements[0]);
      MPI_Get_address(&test_object->target_processor, &displacements[1]);
      MPI_Get_address(&test_object->local, &displacements[2]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      displacements[2] -= test_object_ptr;
      MPI_Datatype types[3] = {dt1.get(), dt2.get(), dt3.get()};
      MPI_Type_create_struct(3, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };

  // Local edge descriptors
  template <typename A, typename B>
  struct make_mpi_datatype<boost::detail::edge_desc_impl<A, B> > 
    : make_mpi_datatype_base 
  {
    make_mpi_datatype<B> dt1;
    make_mpi_datatype<B> dt2;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1(), dt2() {
      int blocklengths[2] = {1, 1};
      MPI_Aint displacements[2];
      char dummy;
      boost::detail::edge_desc_impl<A, B> *test_object
	= (boost::detail::edge_desc_impl<A, B>*)(&dummy);
      MPI_Aint test_object_ptr;
      MPI_Get_address(test_object, &test_object_ptr);
      MPI_Get_address(&test_object->m_source, &displacements[0]);
      MPI_Get_address(&test_object->m_target, &displacements[1]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      MPI_Datatype types[2] = {dt1.get(), dt2.get()};
      MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };

  template <typename Hd>
  struct make_mpi_datatype<boost::tuples::cons<Hd, boost::tuples::cons<boost::no_property, 
								       boost::tuples::null_type> > > 
    : make_mpi_datatype_base {
    make_mpi_datatype<Hd> dt_hd;
    make_mpi_datatype(): dt_hd() {}
    MPI_Datatype get() const {return dt_hd.get();}
  };
}

#endif // BOOST_GRAPH_DISTRIBUTED_ADJACENCY_LIST_HPP
