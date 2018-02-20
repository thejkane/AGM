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

#ifndef AMPLUSPLUS_OBJECT_BASED_ADDRESSING_HPP
#define AMPLUSPLUS_OBJECT_BASED_ADDRESSING_HPP

#include <am++/basic_coalesced_message_type.hpp>
#include <am++/dummy_buffer_sorter.hpp>
#include <am++/detail/factory_wrapper.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range.hpp>
//#include <boost/intrusive/detail/math.hpp>
#include <boost/intrusive/detail/utilities.hpp>
#include <boost/shared_container_iterator.hpp>
#include <boost/format.hpp>
#include <boost/assert.hpp>
#include <iostream>
#include <utility>
#include <stdio.h>

// #define DISABLE_SELF_SEND_CHECK // For debugging

namespace amplusplus {

  namespace detail {
    template <typename Handler>
    class oba_handler {
      public:
#ifndef BOOST_NO_RVALUE_REFERENCES
      oba_handler(const Handler& handler): handler(handler) {}
      oba_handler(Handler&& h = Handler()): handler(std::move(h)) {}
#else
      oba_handler(const Handler& handler = Handler()): handler(handler) {}
#endif
#ifndef BOOST_NO_DEFAULTED_FUNCTIONS
      oba_handler(const oba_handler&) = default;
      oba_handler(oba_handler&&) = default;
      oba_handler& operator=(oba_handler&&) = default;
      oba_handler& operator=(const oba_handler&) = default;
#endif

      template <typename RankType, typename Data>
      void operator()(RankType& /*rank*/, const Data& data) {
        handler(data);
      }

      const Handler& get_handler() const {return handler;}

      private:
      Handler handler;

      template <typename, typename, typename, typename, typename, typename>
      friend class object_based_addressing;
    };

    struct host_based_routing_tag {}; // Used to mark special constructor in OBA that's used for routing
  }

  namespace detail {
    template <typename OBA, typename RankType, typename Handler>
    struct oba_dest_hbr_handler;

    class all_ranks_except_me: public valid_rank_set_base {
      transport::rank_type rank, size;

      public:
      all_ranks_except_me(transport::rank_type rank, transport::rank_type size): rank(rank), size(size) {}
      bool is_valid(transport::rank_type r) const {return r < size && r != rank;}
      transport::rank_type count() const {return size - 1;}
      transport::rank_type rank_from_index(transport::rank_type idx) const {return idx < rank ? idx : idx + 1;}
    };

  }

  template <typename Arg, typename Handler, typename OwnerMap, typename CoalescingLayerGen = basic_coalesced_message_type_gen, typename BufferSorter = DummyBufferSorter<Arg>, typename SentPartMap = boost::typed_identity_property_map<Arg> >
  class object_based_addressing
    : public CoalescingLayerGen::template inner<Arg, detail::oba_handler<Handler>, BufferSorter >::type
  {
    public:
    typedef typename CoalescingLayerGen::template inner<Arg, detail::oba_handler<Handler>, BufferSorter >::type base_type;
    typedef base_type send_base_type;
    typedef transport::rank_type rank_type;
    typedef Handler handler_type;

#ifdef DISABLE_SELF_SEND_CHECK
    object_based_addressing(const CoalescingLayerGen& coalescing_layer_gen_,
                            const transport& trans,
                            const OwnerMap& owner_,
                            const BufferSorter& bufsrter = DummyBufferSorter<Arg>(),
                            const SentPartMap& sent_part_ = SentPartMap())
      : base_type(coalescing_layer_gen_,
                  trans,
                  boost::make_shared<detail::all_ranks>(trans.size()) /* possible_dests */ ,
                  boost::make_shared<detail::all_ranks>(trans.size()) /* possible_sources */,
                  bufsrter
                  /*Buffer sorter goes here */),
        initialized(true), owner(owner_), sent_part(sent_part_), my_rank(trans.rank()), num_ranks(trans.size())
    {
    }
#else
    object_based_addressing(const CoalescingLayerGen& coalescing_layer_gen_,
                            const transport& trans,
                            const OwnerMap& owner_,
                            const BufferSorter& bufsrter = DummyBufferSorter<Arg>(),
                            const SentPartMap& sent_part_ = SentPartMap())
      : base_type(
          coalescing_layer_gen_,
          trans,
          boost::make_shared<detail::all_ranks_except_me>(trans.rank(), trans.size()) /* possible_dests */ ,
          boost::make_shared<detail::all_ranks_except_me>(trans.rank(), trans.size()) /* possible_sources */,
          bufsrter),
        initialized(true), owner(owner_), sent_part(sent_part_), my_rank(trans.rank()), num_ranks(trans.size())
    {
    }
#endif

    object_based_addressing(detail::host_based_routing_tag,
                            const CoalescingLayerGen& coalescing_layer_gen_,
                            const transport& trans_,
                            const OwnerMap& owner_,
                            const valid_rank_set& possible_dests_,
                            const valid_rank_set& possible_sources_,
                            const BufferSorter& bufsrter = DummyBufferSorter<Arg>(),
                            const SentPartMap& sent_part_ = SentPartMap())
      : base_type(coalescing_layer_gen_, trans_, possible_dests_, possible_sources_, bufsrter), initialized(true), owner(owner_), sent_part(sent_part_), my_rank(trans_.rank()), num_ranks(trans_.size())
    {
    }

#ifndef BOOST_NO_DELETED_FUNCTIONS
    object_based_addressing(const object_based_addressing&) = delete;
#endif

#ifndef BOOST_NO_DEFAULTED_FUNCTIONS
#ifndef BOOST_NO_RVALUE_REFERENCES
    object_based_addressing(object_based_addressing&&) = default;
    object_based_addressing& operator=(object_based_addressing&&) = default;
#endif
#endif

#ifndef BOOST_NO_RVALUE_REFERENCES
    void set_handler(handler_type&& handler_) {
      BOOST_ASSERT (initialized);
      base_type::set_handler(detail::oba_handler<Handler>(std::move(handler_)));
    }
#endif

    void set_handler(const handler_type& handler_) {
      BOOST_ASSERT (initialized);
      base_type::set_handler(detail::oba_handler<Handler>(handler_));
    }

    const handler_type& get_handler() {
      BOOST_ASSERT (initialized);
      return base_type::get_handler().get_handler();
    }

    void send(const Arg& arg) {
      BOOST_ASSERT (initialized);
      rank_type dest = get(owner, arg);
      BOOST_ASSERT (dest < num_ranks);
#ifdef DISABLE_SELF_SEND_CHECK
      send_base_type::send(get(sent_part, arg), dest);
#else
      if (dest == my_rank) {
        get_handler()(get(sent_part, arg));
      } else {
        send_base_type::send(get(sent_part, arg), dest);
      }
#endif
    }

    void send_with_tid(const Arg& arg, int tid) {
      BOOST_ASSERT (initialized);
      rank_type dest = get(owner, arg);
      BOOST_ASSERT (dest < num_ranks);
#ifdef DISABLE_SELF_SEND_CHECK
      send_base_type::send_with_tid(get(sent_part, arg), dest, tid);
#else
      if (dest == my_rank) {
        get_handler()(get(sent_part, arg));
      } else {
        send_base_type::send_with_tid(get(sent_part, arg), dest, tid);
      }
#endif
    }

    transport get_transport() const {return base_type::get_transport();}

    private:
    bool initialized;
    OwnerMap owner;
    SentPartMap sent_part;
    rank_type my_rank, num_ranks;

    template <typename OBA2, typename RankType2, typename Handler2>
    friend struct detail::oba_dest_hbr_handler;

    template <typename Arg2, typename Handler2, typename OwnerMap2, typename Routing2, typename CoalescingLayerGen2, typename BufferSorter2>
    friend class object_based_addressing_dest_hbr;
  };

  template <typename Arg, typename Handler, typename Routing, typename NextHop, typename CoalescingLayerGen = basic_coalesced_message_type_gen, typename BufferSorter = DummyBufferSorter<Arg> >
  class object_based_addressing_dest_hbr;

  namespace detail {
    template <typename OBA, typename RankType, typename Handler>
    struct oba_dest_hbr_handler {
      oba_dest_hbr_handler(OBA& oba, Handler handler = Handler())
        : oba(&oba), my_rank(oba.get_transport().rank()), handler(AMPLUSPLUS_MOVE(handler)) {}
      oba_dest_hbr_handler(): oba(0), my_rank(0), handler() {}

      template <typename Arg>
      void operator()(const Arg& arg) const {
        typedef typename OBA::base_type base_type;
        typedef typename OBA::send_base_type send_base_type;
        RankType o = get(oba->base_type::owner, arg);
        if (o == my_rank) {
          // std::cout << "On " << my_rank << ", got " << arg << " which is mine" << std::endl;
          handler(arg);
        } else {
          RankType nx = oba->routing.next_hop(o);
          // std::cout << "On " << my_rank << ", got " << arg << " destined for " << o << ", routing to " << nx << std::endl;
          oba->send_base_type::send(arg, nx);
        }
      }

      OBA* oba;
      RankType my_rank;
      Handler handler;
    };
  }

  template <typename Arg, typename Handler, typename OwnerMap, typename Routing, typename CoalescingLayerGen, typename BufferSorter >
  class object_based_addressing_dest_hbr
    : public object_based_addressing<Arg, detail::oba_dest_hbr_handler<object_based_addressing_dest_hbr<Arg, Handler, OwnerMap, Routing, CoalescingLayerGen, BufferSorter>, transport::rank_type, Handler>, OwnerMap, CoalescingLayerGen, BufferSorter>
  {
    typedef detail::oba_dest_hbr_handler<object_based_addressing_dest_hbr<Arg, Handler, OwnerMap, Routing, CoalescingLayerGen, BufferSorter>, transport::rank_type, Handler> base_handler_type;
    typedef object_based_addressing<Arg, base_handler_type, OwnerMap, CoalescingLayerGen, BufferSorter> base_type;
    typedef typename base_type::send_base_type send_base_type;
    typedef transport::rank_type rank_type;
    typedef Handler handler_type;

    public:
    object_based_addressing_dest_hbr(transport trans_): base_type(trans_), initialized(false), routing() {}

    object_based_addressing_dest_hbr(const std::pair<CoalescingLayerGen, Routing>& p,
                                     const transport& trans_,
                                     const OwnerMap& owner_,
                                     const BufferSorter& bufsrter = DummyBufferSorter<Arg>())
      : base_type(detail::host_based_routing_tag(), p.first, trans_, owner_, p.second.get_possible_dests(), p.second.get_possible_sources(), bufsrter),
        initialized(true),
        routing(p.second)
    {
    }

    void set_handler(const handler_type& handler_) {
      BOOST_ASSERT (initialized);
      base_type::set_handler(base_handler_type(*this, handler_));
    }

#ifndef BOOST_NO_RVALUE_REFERENCES
    void set_handler(handler_type&& handler_) {
      BOOST_ASSERT (initialized);
      base_type::set_handler(base_handler_type(*this, std::move(handler_)));
    }
#endif

    handler_type& get_handler() {
      BOOST_ASSERT (initialized);
      return base_type::get_handler().handler;
    }

    transport get_transport() const {return base_type::get_transport();}

    void send(const Arg& arg) {
      BOOST_ASSERT (initialized);
#ifdef DISABLE_SELF_SEND_CHECK
      transport::rank_type o = get(base_type::owner, arg);
      transport::rank_type nx = this->routing.next_hop(o);
      // std::cout << "On " << my_rank << ", got " << arg << " destined for " << o << ", routing to " << nx << std::endl;
      send_base_type::send(arg, nx);
#else
      base_type::get_handler()(arg);
#endif
    }

    void send_with_tid(const Arg& arg, int tid) {
      BOOST_ASSERT (initialized);
#ifdef DISABLE_SELF_SEND_CHECK
      transport::rank_type o = get(base_type::owner, arg);
      transport::rank_type nx = this->routing.next_hop(o);
      // std::cout << "On " << my_rank << ", got " << arg << " destined for " << o << ", routing to " << nx << std::endl;
      send_base_type::send_with_tid(arg, nx, tid);
#else
      base_type::get_handler()(arg);
#endif
    }

    private:
    bool initialized;
    Routing routing;
    friend struct detail::oba_dest_hbr_handler<object_based_addressing_dest_hbr<Arg, Handler, OwnerMap, Routing, CoalescingLayerGen, BufferSorter>, rank_type, Handler>;
  };

// Routing schemes for object_based_addressing_dest_hbr

struct no_routing {
  transport::rank_type my_rank;
  transport::rank_type sz;
  no_routing(): my_rank(0), sz(0) {}
  no_routing(transport::rank_type my_rank, transport::rank_type sz): my_rank(transport::rank_type(my_rank)), sz(sz) {}
  transport::rank_type next_hop(transport::rank_type i) const {return i;}
  valid_rank_set get_possible_dests() const {
#ifdef DISABLE_SELF_SEND_CHECK
    return boost::make_shared<amplusplus::detail::all_ranks>(sz);
#else
    return boost::make_shared<amplusplus::detail::all_ranks_except_me>(my_rank, sz);
#endif
  }
  valid_rank_set get_possible_sources() const {
#ifdef DISABLE_SELF_SEND_CHECK
    return boost::make_shared<amplusplus::detail::all_ranks>(sz);
#else
    return boost::make_shared<amplusplus::detail::all_ranks_except_me>(my_rank, sz);
#endif
  }
};

struct ring_routing {
  transport::rank_type my_rank;
  transport::rank_type sz;
  ring_routing(): my_rank(0), sz(0) {}
  ring_routing(transport::rank_type my_rank, transport::rank_type sz): my_rank(transport::rank_type(my_rank)), sz(sz) {}
  transport::rank_type next_hop(transport::rank_type i) const {
    return transport::rank_type(i) == my_rank ? i : transport::rank_type(((my_rank + 1) % sz));
  }
  struct ring_possible_dests: valid_rank_set_base {
    transport::rank_type my_rank, sz;
    ring_possible_dests(transport::rank_type my_rank, transport::rank_type sz): my_rank(my_rank), sz(sz) {}
    bool is_valid(transport::rank_type rank) const {return rank == (my_rank + 1) % sz;}
    transport::rank_type count() const {return 1;}
    transport::rank_type rank_from_index(transport::rank_type i) const {BOOST_ASSERT (i == 0); (void)i; return (my_rank + 1) % sz;}
  };
  struct ring_possible_sources: valid_rank_set_base {
    transport::rank_type my_rank, sz;
    ring_possible_sources(transport::rank_type my_rank, transport::rank_type sz): my_rank(my_rank), sz(sz) {}
    bool is_valid(transport::rank_type rank) const {return rank == (my_rank + sz - 1) % sz;}
    transport::rank_type count() const {return 1;}
    transport::rank_type rank_from_index(transport::rank_type i) const {BOOST_ASSERT (i == 0); (void)i; return (my_rank + sz - 1) % sz;}
  };
  valid_rank_set get_possible_dests() const {
    return boost::make_shared<ring_possible_dests>(my_rank, sz);
  }
  valid_rank_set get_possible_sources() const {
    return boost::make_shared<ring_possible_sources>(my_rank, sz);
  }
};

struct hypercube_routing {
  transport::rank_type my_rank;
  transport::rank_type sz;
  transport::rank_type lg_sz;
  hypercube_routing(): my_rank(0), sz(0), lg_sz(0) {}
  hypercube_routing(transport::rank_type my_rank, transport::rank_type sz)
    : my_rank(transport::rank_type(my_rank)), sz(sz), lg_sz(transport::rank_type(boost::intrusive::detail::floor_log2(sz)))
  {
    if ((sz & (sz - 1)) != 0) { // Not a power of 2
      fprintf(stderr, "Trying to use hypercube routing on non-power-of-2 process count %zu\n", sz);
      abort();
    }
  }
  transport::rank_type next_hop(transport::rank_type i) const {
    // Next hop is gotten by changing the lowest bit that differs between
    // my_rank and i to match the value in i, doing nothing if my_rank == i
    transport::rank_type delta = transport::rank_type(i) ^ my_rank;
    return transport::rank_type(my_rank ^ (delta & -delta));
  }
  struct hypercube_possible_dests: valid_rank_set_base {
    transport::rank_type my_rank, sz, lg_sz;
    hypercube_possible_dests(transport::rank_type my_rank, transport::rank_type sz, transport::rank_type lg_sz): my_rank(my_rank), sz(sz), lg_sz(lg_sz) {}
    bool is_valid(transport::rank_type rank) const {
      transport::rank_type delta = my_rank ^ rank;
      return ((delta & (delta - 1)) == 0); // Is a power of 2
    }
    transport::rank_type count() const {return lg_sz;}
    transport::rank_type rank_from_index(transport::rank_type i) const {return my_rank ^ (transport::rank_type(1) << i);}
  };
  valid_rank_set get_possible_dests() const {return boost::make_shared<hypercube_possible_dests>(my_rank, sz, lg_sz);}
  valid_rank_set get_possible_sources() const {
    return this->get_possible_dests(); // Symmetric
  }
};

struct dissemination_routing {
  transport::rank_type my_rank;
  transport::rank_type sz;
  transport::rank_type lg_sz;
  dissemination_routing(): my_rank(0), sz(0), lg_sz(0) {}
  dissemination_routing(transport::rank_type my_rank, transport::rank_type sz)
    : my_rank(transport::rank_type(my_rank)), sz(sz), lg_sz(transport::rank_type(boost::intrusive::detail::ceil_log2(sz))) {}
  transport::rank_type next_hop(transport::rank_type i) const {
    transport::rank_type delta = (sz + transport::rank_type(i) - my_rank) % sz;
    return transport::rank_type((my_rank + (delta & -delta)) % sz);
  }
  struct dissemination_possible_dests: valid_rank_set_base {
    transport::rank_type my_rank, sz, lg_sz;
    dissemination_possible_dests(transport::rank_type my_rank, transport::rank_type sz, transport::rank_type lg_sz): my_rank(my_rank), sz(sz), lg_sz(lg_sz) {}
    bool is_valid(transport::rank_type rank) const {
      transport::rank_type delta = (rank + sz - my_rank) % sz;
      return ((delta & (delta - 1)) == 0); // Is a power of 2
    }
    transport::rank_type count() const {return lg_sz;}
    transport::rank_type rank_from_index(transport::rank_type i) const {return (my_rank + (transport::rank_type(1) << i)) % sz;}
  };
  struct dissemination_possible_sources: valid_rank_set_base {
    transport::rank_type my_rank, sz, lg_sz;
    dissemination_possible_sources(transport::rank_type my_rank, transport::rank_type sz, transport::rank_type lg_sz): my_rank(my_rank), sz(sz), lg_sz(lg_sz) {}
    bool is_valid(transport::rank_type rank) const {
      transport::rank_type delta = (my_rank + sz - rank) % sz;
      return ((delta & (delta - 1)) == 0); // Is a power of 2
    }
    transport::rank_type count() const {return lg_sz;}
    transport::rank_type rank_from_index(transport::rank_type i) const {return (my_rank + sz - (transport::rank_type(1) << i)) % sz;}
  };
  valid_rank_set get_possible_dests() const {return boost::make_shared<dissemination_possible_dests>(my_rank, sz, lg_sz);}
  valid_rank_set get_possible_sources() const {return boost::make_shared<dissemination_possible_sources>(my_rank, sz, lg_sz);}
};

struct rook_routing {
  transport::rank_type my_rank;
  transport::rank_type side1, side2, sz;
  rook_routing(): my_rank(0), side1(0), side2(0), sz(0) {}
  rook_routing(transport::rank_type my_rank, transport::rank_type sz)
    : my_rank(transport::rank_type(my_rank)), sz(sz)
  {
    transport::rank_type factor = int(sqrt(double(sz)));
    while (sz % factor != 0) --factor; // Must finish at least by when factor == 1
    side1 = factor;
    side2 = sz / factor;
    if (side2 > 3 * side1 && my_rank == 0) {
      std::clog << (boost::format("rook_routing: rectangle sides %d and %d are unbalanced\n") % side1 % side2).str() << std::flush;
    }
  }
  transport::rank_type next_hop(transport::rank_type i) const {
    BOOST_ASSERT (sz != 0);
    if (i / side1 == my_rank / side1) {
      return i;
    } else {
      return transport::rank_type((i / side1) * side1 + my_rank % side1);
    }
  }
  struct rook_possible_dests: valid_rank_set_base {
    transport::rank_type my_rank, side1, side2;
    rook_possible_dests(transport::rank_type my_rank, transport::rank_type side1, transport::rank_type side2): my_rank(my_rank), side1(side1), side2(side2) {}
    bool is_valid(transport::rank_type rank) const {
      return ((rank / side1 == my_rank / side1 || rank % side1 == my_rank % side1) && rank != my_rank);
    }
    transport::rank_type count() const {return side1 + side2 - 2;}
    transport::rank_type rank_from_index(transport::rank_type i) const {
      if (i < side1 - 1) {
        return transport::rank_type((my_rank / side1) * side1 + i + (i >= my_rank % side1));
      } else {
        i -= side1 - 1;
        return transport::rank_type((i + (i >= my_rank / side1)) * side1 + my_rank % side1);
      }
    }
  };
  valid_rank_set get_possible_dests() const {
    BOOST_ASSERT (sz != 0);
    return boost::make_shared<rook_possible_dests>(my_rank, side1, side2);
  }
  valid_rank_set get_possible_sources() const {
    return this->get_possible_dests(); // Symmetric
  }
};

}

#endif // AMPLUSPLUS_OBJECT_BASED_ADDRESSING_HPP
