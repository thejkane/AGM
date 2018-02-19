// Copyright 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_PARALLEL_BFS_HPP
#define BOOST_GRAPH_PARALLEL_BFS_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/detail/thread_support.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/distributed/queue.hpp>
#include <boost/parallel/append_buffer.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/parallel/properties.hpp>
#include <boost/graph/parallel/container_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/graph_utility.hpp> // for map_of_project1st

namespace boost {
  namespace detail {
    /** @brief A unary predicate that decides when to push into a
     *         breadth-first search queue.
     *
     *  This predicate stores a color map that is used to determine
     *  when to push. If it is provided with a key for which the color
     *  is white, it darkens the color to gray and returns true (so
     *  that the value will be pushed appropriately); if the color is
     *  not white, it returns false so that the vertex will be
     *  ignored.
     */
    template<typename ColorMap, typename OwnerMap, typename Rank>
    struct darken_and_push
    {
      typedef typename property_traits<ColorMap>::key_type argument_type;
      typedef bool result_type;

      explicit darken_and_push(const ColorMap& color, const OwnerMap& owner, Rank rank)
        : color(color), owner(owner), rank(rank) { }

      bool operator()(const argument_type& value) const
      {
        typedef color_traits<typename property_traits<ColorMap>::value_type> Color;

        if (get(color, value) == Color::white()) {
          if (get(owner, value) == rank)
            return exchange(color, value, Color::white(), Color::gray());
          else
            return true; // Could also optionally perform a cache() here if we had thread-safe ghost cells 
        }
        else 
          return false;
      }

    private:
      ColorMap color;
      OwnerMap owner;
      Rank     rank;
    };

    template<typename IndexMap>
    struct has_not_been_seen
    {
      typedef bool result_type;

      has_not_been_seen() { }

      has_not_been_seen(std::size_t n, IndexMap index_map)
        : seen(n), index_map(index_map) {}

      template<typename Key>
      result_type operator()(Key key)
      {
        bool result = seen[get(index_map, key)];
        seen[get(index_map, key)] = true;
        return !result;
      }

      void swap(has_not_been_seen& other)
      {
        using std::swap;
        swap(seen, other.seen);
        swap(index_map, other.index_map);
      }

    private:
      dynamic_bitset<> seen;
      IndexMap index_map;
    };

    template<typename IndexMap>
    inline void
    swap(has_not_been_seen<IndexMap>& x, has_not_been_seen<IndexMap>& y)
    {
      x.swap(y);
    }

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
    struct hit_rate_counter {
      amplusplus::detail::atomic<unsigned long> hits, tests, local_msgs, remote_msgs;
      hit_rate_counter(): hits(0), tests(0), local_msgs(0), remote_msgs(0) {}
      void clear() {hits.store(0); tests.store(0);}
      void test() {tests.fetch_add(1);}
      void hit() {hits.fetch_add(1);}
      void local_msg() {local_msgs.fetch_add(1);}
      void remote_msg() {remote_msgs.fetch_add(1);}
    };
#else
    struct hit_rate_counter {
      void clear() {}
      void test() {}
      void hit() {}
    };
#endif

  } // end namespace detail

  namespace graph {
    namespace distributed {

    template <class Graph, class BFSVisitor = bfs_visitor<>, 
              class ColorMap = 
                two_bit_color_map<typename property_map<Graph, vertex_index_t>::const_type>,
              class Buffer =    
                distributed_queue<typename property_map<Graph, vertex_owner_t>::const_type,
                                  append_buffer<typename graph_traits<Graph>::vertex_descriptor>,
                                  boost::detail::darken_and_push<ColorMap, 
                                                                 typename property_map<Graph, vertex_owner_t>::const_type,
                                                                 typename Graph::rank_type> > >
    class breadth_first_search
    {
      typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename property_traits<ColorMap>::value_type ColorValue;
      typedef color_traits<ColorValue> Color;
      typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap; 
      typedef typename Graph::rank_type rank_type;

    public:

      breadth_first_search(Graph& g) 
        : transport(g.transport()), g(g), vis(make_bfs_visitor(null_visitor())), 
          color(num_vertices(g), get(vertex_index, g)),
          Q(new Buffer(g.transport(), 
                       get(vertex_owner, g), 
                       make_shared<append_buffer<Vertex> >(),
                       make_shared<append_buffer<Vertex> >(),
                       boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), transport.rank())))
      { }

      breadth_first_search(Graph& g, BFSVisitor vis) 
        : transport(g.transport()), g(g), vis(vis),
          color(num_vertices(g), get(vertex_index, g)),
          Q(new Buffer(g.transport(), 
                       get(vertex_owner, g), 
                       make_shared<append_buffer<Vertex> >(),
                       make_shared<append_buffer<Vertex> >(),
                       boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), transport.rank())))
      { }

      breadth_first_search(Graph& g, BFSVisitor vis, ColorMap color) 
        : transport(g.transport()), g(g), vis(vis), color(color),
          Q(new Buffer(g.transport(), 
                       get(vertex_owner, g), 
                       make_shared<append_buffer<Vertex> >(),
                       make_shared<append_buffer<Vertex> >(),
                       boost::detail::darken_and_push<ColorMap, OwnerMap, rank_type>(color, get(vertex_owner, g), transport.rank())))
      { }

      breadth_first_search(Graph& g, BFSVisitor vis, ColorMap color,
                           shared_ptr<Buffer> Q) 
        : transport(g.transport()), g(g), vis(vis), color(color), Q(Q) 
      { }

      void set_source(Vertex s) { source = s; } // for threaded execution
      void operator() (int tid) { run_no_init(source, tid); }

      void run(Vertex s, int tid = 0) 
      { 
        if (tid == 0) {
          set_property_map_role(vertex_color, color);
          color.set_consistency_model(0);

          BGL_FORALL_VERTICES_T(v, g, Graph) {
            vis.initialize_vertex(v, g);
            put(color, v, Color::white());
          }
        }

        run_no_init(s, tid);
      }

      void run_no_init(Vertex s, int tid = 0) 
      { 
        AMPLUSPLUS_WITH_THREAD_ID(tid) {

          if (tid == 0)
            t_bar.reset(new amplusplus::detail::barrier(transport.get_nthreads()));
          
          // This scoped_epoch acts as a temporary barrier until we can be sure t_bar is initialized 
          {
            amplusplus::scoped_epoch epoch(transport);

            if (s != graph_traits<Graph>::null_vertex() && // Lets us push a source value we don't intend to use
                get(vertex_owner, g, s) == transport.rank() && tid == 0) {
              Q->push(s);
            }
          }

	  {
	    // This epoch may not be necessary, but we put the whole thing in an epoch just in case the put below will want to send some messages (it should not since we ensure we are the owner of the vertex, but I am not sure if we currently forbid sending messages from puts that are local).
	    amplusplus::scoped_epoch epoch(transport);
	    if (s != graph_traits<Graph>::null_vertex() && // Lets us push a source value we don't intend to use
		get(vertex_owner, g, s) == transport.rank() && tid == 0) {
	    // We cannot update the color in the same epoch where we push the source onto the queue.  If we are using the darken_and_push predicate for the queue, it may prevent the push_handler of the queue from actually putting source on the queue if we prematurly set the color to black.
	      put(color, s, Color::black()); 
	    }
	  }

          unsigned long all_queues_empty;

          while (true) {
            t_bar->wait();
            if (tid == 0) Q->swap();
            t_bar->wait();

            const unsigned long my_queue_empty = (tid == 0 && Q->local_empty()) ? 1 : 0;
            {
              amplusplus::scoped_epoch_value epoch(transport, my_queue_empty, all_queues_empty);
            
              typename Buffer::size_type nthr = transport.get_nthreads();
              for (typename Buffer::size_type i = tid ; i < Q->size() ; i += nthr) {
                assert(i < Q->size());
                Vertex v = (*Q)[i];

                BGL_FORALL_ADJ_T(v, u, g, Graph) { Q->push(u); }
              
                put(color, v, Color::black());
              }
            }
            if (all_queues_empty == transport.size()) break;

            t_bar->wait();
            if (tid == 0) Q->clear();
          }
        }
      }

    private:
      amplusplus::transport&                  transport;
      const Graph&                            g;
      BFSVisitor                              vis;
      ColorMap                                color;
      shared_ptr<Buffer>                      Q;
      shared_ptr<amplusplus::detail::barrier> t_bar;
      Vertex                                  source; // for threaded execution
    };

    //
    // Fully-asynchronous, message-driven BFS 
    //
    template <class Graph, class DistanceMap,
              class MessageGenerator = 
                amplusplus::simple_generator<amplusplus::basic_coalesced_message_type_gen> >
    class async_breadth_first_search {
      
      typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename property_traits<DistanceMap>::value_type Distance;

      typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap;
      typedef boost::parallel::map_of_project1st<OwnerMap> MsgOwnerMap;

      typedef std::pair<Vertex, Distance> msg_value_type;

    public:

      async_breadth_first_search(Graph& g, DistanceMap& distance,
				 MessageGenerator message_gen = 
		                   MessageGenerator(amplusplus::basic_coalesced_message_type_gen(1 << 12)))
	: dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_value_type>(), 0)),
	  transport(g.transport()), g(g), distance(distance), owner(get(vertex_owner, g)), msg_owner(owner),
	  msg(message_gen, transport, msg_owner, 
	      amplusplus::idempotent_combination(boost::parallel::minimum<Distance>(), 
						 std::numeric_limits<Distance>::max()))
      {
	msg.set_handler(msg_handler(*this));
      }

      void run(Vertex v) { 
        AMPLUSPLUS_WITH_THREAD_ID(0) {

	  amplusplus::scoped_epoch epoch(transport);
	  if (get(owner, v) == transport.rank())
	    msg.send(msg_value_type(v, 0)); 
	}
      }

      std::pair<unsigned long, unsigned long>
      get_cache_stats() {
	boost::parallel::all_reduce<unsigned long, std::plus<unsigned long> > 
	  r(transport, std::plus<unsigned long>());
	
	unsigned long total_tests = r(msg.get()->counters.tests.load());
	unsigned long total_hits = r(msg.get()->counters.hits.load());
	
	return std::make_pair(total_tests, total_hits);
      }

#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      std::pair<unsigned long, unsigned long>
      get_msg_stats() {
        boost::parallel::all_reduce<unsigned long, std::plus<unsigned long> > 
          r(transport, std::plus<unsigned long>());
        
        unsigned long total_tests = r(counters.tests.load());
        unsigned long total_hits = r(counters.hits.load());
        
//      boost::parallel::all_reduce<size_t, std::plus<size_t> >
//        r2(transport, std::plus<size_t>());
//      for (size_t i = 0 ; i < distance_updates.size() ; ++i) {
//        size_t t = r2(distance_updates[i]);
//        if (transport.rank() == 0) std::cout << i << ": " << t << std::endl; 
//        if (t == 0) break;
//      }

        unsigned long total_local_msgs = r(counters.local_msgs.load());
        unsigned long total_remote_msgs = r(counters.remote_msgs.load());
        
        if (transport.rank() == 0)
          std::cout << "  " << total_local_msgs << " local messages and " << total_remote_msgs << " remote messages\n";
  
        return std::make_pair(total_tests, total_hits);
      }
#endif

    private:

      struct msg_handler;

      typedef typename MessageGenerator::template call_result<msg_value_type, msg_handler, MsgOwnerMap,
							      amplusplus::idempotent_combination_t<boost::parallel::minimum<Distance>, 
				                                                Distance> >::type
        message_type;

      const int dummy_first_member_for_init_order; // Unused
      amplusplus::transport& transport;
      const Graph& g;
      DistanceMap distance;
      OwnerMap owner;
      MsgOwnerMap msg_owner;

      message_type msg;
      boost::detail::hit_rate_counter counters;
//       std::vector<size_t> distance_updates;
    };


    template <typename Graph, typename DistanceMap, typename MessageGenerator>
    struct async_breadth_first_search<Graph, DistanceMap, MessageGenerator>::msg_handler {
      
      msg_handler() : self(NULL) {}
      msg_handler(async_breadth_first_search& self) : self(&self) {}
      
      void operator() (const msg_value_type& data) const
      {
        Vertex v = data.first;
        Distance old_dist;
        Distance new_dist = data.second;

        do {
          old_dist = get(self->distance, v);
          self->counters.test();
          if (new_dist >= old_dist) return;
        } while (!exchange(self->distance, v, old_dist, new_dist));

        self->counters.hit();
//      self->distance_updates[new_dist]++;

        // TODO: Explore all out edges of v
        BGL_FORALL_ADJ_T(v, u, self->g, Graph) {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
          if (get(self->owner, u) == self->transport.rank())
            self->counters.local_msg();
          else
            self->counters.remote_msg();
#endif
 
	  self->msg.send(std::make_pair(u, new_dist + 1));
	}
      }

    protected:
      async_breadth_first_search* self;
    };


    template <class Graph, class VisitedMap,
              class MessageGenerator = 
                amplusplus::simple_generator<amplusplus::basic_coalesced_message_type_gen> >
    class parallel_search {
      
      typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename property_traits<VisitedMap>::value_type VisitedValue;

      typedef typename property_map<Graph, vertex_owner_t>::const_type OwnerMap;

    public:

      parallel_search(Graph& g, const VisitedMap& visited, VisitedValue visited_value,
		      MessageGenerator message_gen = 
		      MessageGenerator(amplusplus::basic_coalesced_message_type_gen(1 << 12)))
	: dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<Vertex>(), 0)),
	  transport(g.transport()), g(g), visited(visited), owner(get(vertex_owner, g)),
	  visited_value(visited_value),
	  msg(message_gen, transport, owner, amplusplus::duplicate_removal(boost::parallel::identity<Vertex>()))
      {
	msg.set_handler(msg_handler(*this));
      }

      void run(Vertex v) { 
        AMPLUSPLUS_WITH_THREAD_ID(0) {

	  amplusplus::scoped_epoch epoch(transport);
	  if (v != graph_traits<Graph>::null_vertex() && get(owner, v) == transport.rank())
	    msg.send(v); 
	}
      }

      std::pair<unsigned long, unsigned long>
      get_cache_stats() {
	boost::parallel::all_reduce<unsigned long, std::plus<unsigned long> > 
	  r(transport, std::plus<unsigned long>());
	
	unsigned long total_tests = r(msg.get()->counters.tests.load());
	unsigned long total_hits = r(msg.get()->counters.hits.load());
	
	return std::make_pair(total_tests, total_hits);
      }

    private:

      struct msg_handler;

      typedef typename MessageGenerator::template call_result<Vertex, msg_handler, OwnerMap,
				     amplusplus::duplicate_removal_t<boost::parallel::identity<Vertex> > >::type
        message_type;

      const int dummy_first_member_for_init_order; // Unused
      amplusplus::transport& transport;
      const Graph&           g;
      const VisitedMap&      visited;
      OwnerMap               owner;
      VisitedValue           visited_value;

      message_type msg;
    };


    template <typename Graph, typename DistanceMap, typename MessageGenerator>
    struct parallel_search<Graph, DistanceMap, MessageGenerator>::msg_handler {
      
      msg_handler() : self(NULL) {}
      msg_handler(parallel_search& self) : self(&self) {}
      
      void operator() (const Vertex& v) const
      {
        if (get(self->visited, v) == self->visited_value) return; 

        // Could use atomic version here but it doesn't matter if multiple messages proceed
        put(self->visited, v, self->visited_value); 

	BGL_FORALL_ADJ_T(v, u, self->g, Graph) {
	  self->msg.send(u);
	}
      }

    protected:
      parallel_search* self;
    };
  

} } } // end namespace boost::graph::distributed

#endif // BOOST_GRAPH_PARALLEL_BFS_HPP
