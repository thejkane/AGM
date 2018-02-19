// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Jeremiah Willcock
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_PERMUTE_GRAPH_HPP
#define BOOST_GRAPH_PERMUTE_GRAPH_HPP

#include <boost/mpi.hpp>
#include <boost/ref.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <vector>
#include <am++/am++.hpp>
#include <am++/counter_coalesced_message_type.hpp>
#include <am++/message_type_generators.hpp>
#include <stdint.h>

namespace boost {
  namespace detail {
    template <typename Out>
    class iter_write_handler {
      Out* out;
      iter_write_handler(): out(0) {}
      explicit iter_write_handler(Out& out): out(&out) {}
      template <typename T>
      void operator()(const T& val) const {
	*(*out)++ = val;
      }
    };

    template <typename OutVec>
    class vector_write_handler {
      //    private:
      //      mutable boost::mutex mtx_;
    public:
      OutVec* out;
      vector_write_handler() : out(0) {}
      vector_write_handler(const vector_write_handler& other) : out(other.out) {}
      explicit vector_write_handler(OutVec* pout): out(pout) {}

      vector_write_handler& operator=(const vector_write_handler& other) {
	if (this != &other) {
	  this->out = other.out;
	}
	return *this;
      }

      template <typename T>
      void operator()(const T& val) const {
	const int tid = amplusplus::detail::get_thread_id();
	out[tid].push_back(val);
      }
    };


    template <typename Dist>
    class distrib_to_pair_owner_map_t {
      Dist dist;
      public:
      typedef typename Dist::size_type dist_input_type;
      //      typedef std::pair<dist_input_type, dist_input_type> key_type;
      typedef typename Dist::rank_type value_type;

      explicit distrib_to_pair_owner_map_t(const Dist& dist): dist(dist) {}

      template <typename key_type>
      friend value_type
      get(const distrib_to_pair_owner_map_t& o, const key_type& v) {
        return o.dist(v.first);
      }
    };

    template <typename Dist>
    distrib_to_pair_owner_map_t<Dist> distrib_to_pair_owner_map(const Dist& d) {
      return distrib_to_pair_owner_map_t<Dist>(d);
    }

    struct dummy_flip {template <typename T> T operator()(const T&) const {abort();}};
  }

  // TODO: Add message generator here
  template <typename InIter, typename Flip, typename OutIter, typename Distribution,
            typename MessageGenerator>
  void distribute(InIter begin, InIter end, const Flip& flip, OutIter out,
                  const Distribution& dist,
                  amplusplus::transport& trans,
                  MessageGenerator msg_gen) 
  {
    typedef typename std::iterator_traits<InIter>::value_type value_type;
    typedef detail::iter_write_handler<OutIter> iter_write_handler;

    amplusplus::register_mpi_datatype<value_type>();
    
    typedef typename MessageGenerator::template call_result<value_type, iter_write_handler,
      detail::distrib_to_pair_owner_map_t<Distribution>, 
                                                            amplusplus::no_reduction_t>::type 
      write_msg_type;
    write_msg_type write_msg(msg_gen, trans, detail::distrib_to_pair_owner_map<Distribution>(dist), 
			     amplusplus::no_reduction);

    write_msg.set_handler(iter_write_handler(out));
    {
      amplusplus::scoped_epoch ep(trans);
      for (; begin != end; ++begin) {
        if (begin->first == (begin->second).first) continue;
        write_msg.send(*begin);
        if (!boost::is_same<Flip, detail::dummy_flip>::value) {
          write_msg.send(flip(*begin));
        }
      }
    }
  }

  template <typename InIter, typename Flip, typename SelfLoop, typename MessageType>
  void concurrent_distribute(InIter begin, InIter end, const Flip& flip,
			     const SelfLoop is_selfloop,
                  amplusplus::transport& trans,
                  MessageType& write_msg) {
    amplusplus::scoped_epoch ep(trans);
    for (; begin != end; ++begin) {

      if (is_selfloop(*begin))
	continue;

      write_msg.send(*begin);
      if (!boost::is_same<Flip, detail::dummy_flip>::value) {
	write_msg.send(flip(*begin));
      }
    }
  }


  // TODO: Overload for default flip with MessageGenerator specified?

  // Default MessageGenerator and Flip
  template <typename InIter, typename OutIter, typename Distribution>
  void distribute(InIter begin, InIter end, OutIter out,
                  const Distribution& dist,
                  const amplusplus::transport& trans) 
  {
    typedef typename amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, 
                                                   amplusplus::hypercube_routing>
      MessageGenerator;

    distribute(begin, end, detail::dummy_flip(), out, dist, trans,
               MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12), 
                                amplusplus::hypercube_routing(trans.rank(), trans.size()))); 
  }

  // Default MessageGenerator
  template <typename InIter, typename Flip, typename OutIter, typename Distribution>
  void distribute(InIter begin, InIter end, const Flip& flip, OutIter out,
                  const Distribution& dist,
                  const amplusplus::transport& trans) 
  {
    typedef typename amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, 
                                                   amplusplus::hypercube_routing>
      MessageGenerator;

    distribute(begin, end, flip, out, dist, trans,
               MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12), 
                                amplusplus::hypercube_routing(trans.rank(), trans.size()))); 
  }

  // Default MessageGenerator for distribute below 
  template <typename ValueType, typename Flip, typename Gen, typename OutIter, 
            typename Distribution>
  void distribute(Gen data_gen, uintmax_t start, uintmax_t end,
                  const Flip& flip,
                  OutIter out,
                  const Distribution& dist,
                  const amplusplus::transport& trans)
  {
    typedef typename amplusplus::routing_generator<amplusplus::counter_coalesced_message_type_gen, 
                                                   amplusplus::rook_routing>
      MessageGenerator;

    distribute<ValueType, Flip, Gen, OutIter, Distribution, MessageGenerator>
      (data_gen, start, end, flip, out, dist, trans, 
       MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12), 
                        amplusplus::rook_routing(trans.rank(), trans.size()))); 
  }

  template <typename ValueType, typename Flip, typename Gen, typename OutIter, 
            typename Distribution, typename MessageGenerator> 
  void distribute(Gen data_gen, uintmax_t start, uintmax_t end,
                  const Flip& flip,
                  OutIter out,
                  const Distribution& dist,
                  const amplusplus::transport& trans,
                  MessageGenerator msg_gen) 
  {
    typedef ValueType value_type;

    typedef detail::iter_write_handler<OutIter> iter_write_handler;

    amplusplus::register_mpi_datatype<value_type>();
    
    typedef typename MessageGenerator::template call_result<value_type, iter_write_handler, 
							    detail::distrib_to_pair_owner_map_t<Distribution>, 
							    amplusplus::no_reduction_t>::type 
      write_msg_type;
    write_msg_type write_msg(msg_gen, trans, detail::distrib_to_pair_owner_map(dist), amplusplus::no_reduction);

    write_msg.set_handler(iter_write_handler(out));
    boost::scoped_array<value_type> buf(new value_type[1 << 16]);
    {
      amplusplus::scoped_epoch ep(trans);
      for (; start < end; start += (1 << 16)) {
        uintmax_t local_end = (std::min)(end, start + (1 << 16));
        data_gen(start, local_end, buf.get());
        for (int i = 0; i < int(local_end - start); ++i) {
          write_msg.send(buf[i]);
          if (!boost::is_same<Flip, detail::dummy_flip>::value) {
            write_msg.send(flip(buf[i]));
          }
        }
      }
    }
  }

  template <typename ValueType, typename Gen, typename OutIter, typename Distribution>
  void distribute(Gen data_gen, uintmax_t start, uintmax_t end,
                  OutIter out,
                  const Distribution& dist,
                  const amplusplus::transport& trans) {
    distribute<ValueType>(data_gen, start, end, out, detail::dummy_flip(), dist, trans);
  }

  template <typename Graph, typename RandomGenerator, typename Distribution>
  void rand_dist(typename graph_traits<Graph>::vertices_size_type nverts, 
                 RandomGenerator& gen_orig, Distribution dist, 
                 const boost::mpi::communicator& comm, 
                 std::vector<typename graph_traits<Graph>::vertices_size_type>& result) {

    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;

    boost::mpi::communicator local_comm(comm, boost::mpi::comm_duplicate);
    int rank = local_comm.rank(), size = local_comm.size();
    typename RandomGenerator::split_iterator spliti = gen_orig.split_off_n(size).first;
    std::advance(spliti, rank);
    RandomGenerator gen = *spliti;
    vertices_size_type chunk_size = (nverts + size - 1) / size;
    vertices_size_type my_start = rank * chunk_size;
    vertices_size_type my_end = (std::min)(nverts, (rank + 1) * chunk_size);
    std::vector<std::vector<vertices_size_type> > outbuf(size), result_temp(size);
    boost::uniform_int<> ui(0, size - 1);
    for (vertices_size_type i = my_start; i < my_end; ++i) {
      int dest = ui(gen);
      outbuf[dest].push_back(i);
    }
    boost::mpi::all_to_all(local_comm, outbuf, result_temp);
    result.clear();
    for (int i = 0; i < size; ++i) {
      result.insert(result.end(), result_temp[i].begin(), result_temp[i].end());
    }
    vertices_size_type my_pos_in_final_result = boost::mpi::scan(local_comm, result.size(), std::plus<size_t>()) - result.size();
    for (vertices_size_type i = result.size() - 1; i > 0; --i) {
      vertices_size_type swap_pos = boost::uniform_int<vertices_size_type>(0, i)(gen);
      if (swap_pos != i) {
        std::swap(result[swap_pos], result[i]);
      }
    }
    for (vertices_size_type i = 0; i < size; ++i) {outbuf[i].clear(); result_temp[i].clear();}
    for (vertices_size_type i = 0; i < result.size(); ++i) {
      outbuf[dist(my_pos_in_final_result + i)].push_back(result[i]);
    }
    boost::mpi::all_to_all(local_comm, outbuf, result_temp);
    result.clear();
    for (vertices_size_type i = 0; i < size; ++i) {
      result.insert(result.end(), result_temp[i].begin(), result_temp[i].end());
    }
  }

  template <bool CopyAndFlip, typename Graph, typename InputEdgeIterator, typename Distribution, 
            typename PermutationMap, typename OutputEdgeIterator, typename Transport>
  void permute_and_redistribute_ring(InputEdgeIterator begin, InputEdgeIterator end, 
                                     Distribution dist, PermutationMap permutation_map, 
                                     OutputEdgeIterator out, Transport& trans,
                                     typename graph_traits<Graph>::edges_size_type chunk_size = 1 << 17) {

    typedef typename graph_traits<Graph>::edges_size_type edges_size_type;

    MPI_Comm comm;
    MPI_Comm_dup(trans.template downcast_to_impl<amplusplus::mpi_transport_event_driven>()->get_mpi_communicator(), &comm);

#ifdef BOOST_MPI
    boost::mpi::communicator mpi_comm(comm, boost::mpi::comm_attach);
#endif

    unsigned int rank = trans.rank(), size = trans.size();
    typedef typename std::iterator_traits<InputEdgeIterator>::value_type edge_type;
    int next = (rank + 1) % size, prev = (rank + size - 1) % size;
    // Bits in statuses (from LSB):
    // 3: Wrote flipped final edge (CopyAndFlip only)
    // 2: Wrote final edge
    // 1: Permuted target
    // 0: Permuted source
    std::vector<char> statuses;
    std::vector<edge_type> sendbuf;
    int round = 0;
    edges_size_type edges_output = 0;
    edges_size_type edges_generated = 0;
    while (true) {
      edges_size_type sendcount = statuses.size(), recvcount;

      // TODO: Switch to isend/irecv to use entirely Boost.MPI calls
      MPI_Sendrecv(&sendcount, 1, MPI_UNSIGNED_LONG, next, 0,
                   &recvcount, 1, MPI_UNSIGNED_LONG, prev, 0,
                   comm, MPI_STATUS_IGNORE);
      std::vector<char> instatuses(recvcount);
      std::vector<edge_type> recvbuf(recvcount);
#ifdef BOOST_MPI
      boost::mpi::request reqs[4];
      reqs[0] = mpi_comm.isend(next, 1, statuses);
      reqs[1] = mpi_comm.isend(next, 2, sendbuf);
      reqs[2] = mpi_comm.irecv(prev, 1, instatuses);
      reqs[3] = mpi_comm.irecv(prev, 2, recvbuf);
#else
      MPI_Request reqs[4];
      MPI_Isend(&statuses[0], statuses.size(), MPI_CHAR, next, 1, comm, &reqs[0]);
      MPI_Isend(&sendbuf[0], sendbuf.size() * sizeof(edge_type), MPI_BYTE, next, 2, comm, &reqs[1]);
      MPI_Irecv(&instatuses[0], instatuses.size(), MPI_CHAR, prev, 1, comm, &reqs[2]);
      MPI_Irecv(&recvbuf[0], recvbuf.size() * sizeof(edge_type), MPI_BYTE, prev, 2, comm, &reqs[3]);
#endif

      ++round;
      std::vector<edge_type> generated_edges;
      // Generate edges if there are any left
      if (begin != end) {
        for (edges_size_type count = 0;
             count < (chunk_size) / (CopyAndFlip ? 2 : 1) && begin != end;
             ++count, ++begin) {
          edge_type e = *begin;
          generated_edges.push_back(e);
          ++edges_generated;
        }
      }

#ifdef BOOST_MPI
      boost::mpi::wait_all(&reqs[0], &reqs[4]);
#else
      MPI_Waitall(4, &reqs[0], MPI_STATUS_IGNORE);
#endif

      statuses.clear();
      sendbuf.clear();
      // Handle incoming edges
      for (edges_size_type i = 0; i < recvcount + generated_edges.size(); ++i) {
        char stat = i < recvcount ? instatuses[i] : 0;
        edge_type e = i < recvcount ? recvbuf[i] : generated_edges[i - recvcount];
        if (!(stat & 1) && dist(e.first) == rank) {
          // Source not permuted and I own it
          e.first = get(permutation_map, e.first);
          stat |= 1;
        }
        if (!(stat & 2) && dist(e.second) == rank) {
          // Target not permuted and I own it
          e.second = get(permutation_map, e.second);
          stat |= 2;
        }
        if ((stat & 7) == 3 && dist(e.first) == rank) {
          // Edge fully permuted, not yet written unflipped, and I own it
          *out++ = e;
          stat |= 4;
          ++edges_output;
        }
        if (CopyAndFlip && (stat & 11) == 3 && dist(e.second) == rank) {
          // Edge fully permuted, not yet written flipped, and I own it
          std::swap(e.first, e.second);
          *out++ = e;
          stat |= 8;
          ++edges_output;
        }
        if (stat != (CopyAndFlip ? 15 : 7)) {
          // Edge not fully processed and saved yet
          // Send it to the next process
          statuses.push_back(stat);
          sendbuf.push_back(e);
        }
      }

      boost::parallel::all_reduce<edges_size_type, std::plus<edges_size_type> > r(trans, std::plus<edges_size_type>());
      if( r(sendbuf.size() + edges_generated) == 0)
        break;

      edges_generated = 0;
    }
  }

  template <typename Distribution, typename T>
  struct local_map {
    boost::shared_ptr<std::vector<T> > data;
    const Distribution& dist;
    local_map(boost::shared_ptr<std::vector<T> > data, 
              const Distribution& dist)
      : data(data), dist(dist) {}
  };

  template <typename Distribution, typename T>
  inline size_t get(const local_map<Distribution, T>& lm, T idx) {
    T local_idx = lm.dist.local(idx);
    assert (local_idx < lm.data->size());
    return (*lm.data)[local_idx];
  }

  // Assume InputEdgeIterator::value_type is std::pair<vertex_index>
  // RandomGenerator must be splittable
  template <bool CopyAndFlip, typename Graph, typename InputEdgeIterator, 
            typename RandomGenerator, typename Distribution, typename OutputEdgeIterator>
  void permute_and_redistribute_simple(typename graph_traits<Graph>::vertices_size_type nverts, 
                                       InputEdgeIterator begin, InputEdgeIterator end, 
                                       RandomGenerator& gen, const Distribution& dist, 
                                       OutputEdgeIterator out, 
                                       const boost::mpi::communicator& mpi_comm) {

    // TODO: Static assert graph_traits<Graph>::vertices_size_type is convertible to 
    //       iterator_traits<InputEdgeIterator>::value_type::first_type

    typedef typename std::iterator_traits<InputEdgeIterator>::value_type::first_type vertex_index;
    typedef typename graph_traits<Graph>::edges_size_type edges_size_type;
    
    //
    // Need a vertex_index map to build a property map
    //

    boost::shared_ptr<std::vector<vertex_index> > perm(new std::vector<vertex_index>());
    rand_dist<Graph>(nverts, gen, dist, mpi_comm, *perm);
    local_map<Distribution, edges_size_type> lm(perm, dist);
    permute_and_redistribute_ring<CopyAndFlip, Graph>(begin, end, dist, lm, out, mpi_comm);
  }
}

#endif // BOOST_GRAPH_PERMUTE_GRAPH_HPP
