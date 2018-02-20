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

#ifndef AMPLUSPLUS_MESSAGE_TYPE_GENERATORS_HPP
#define AMPLUSPLUS_MESSAGE_TYPE_GENERATORS_HPP

#include <boost/parameter.hpp>
#include <boost/assert.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/typeof/typeof.hpp>
#include <am++/transport.hpp>
#include <am++/object_based_addressing.hpp>
#include <am++/reductions.hpp>
#include <am++/dummy_buffer_sorter.hpp>
#include <am++/detail/typed_in_place_factory_owning.hpp>


namespace amplusplus {

  // Generator concept
  // concept MessageTypeGenerator<typename Gen> {
  //   template <typename Arg, typename Handler>
  //   struct Gen::type_normal {
  //     typename type;
  //     requires MessageType<type>;
  //   };
  //   template <typename Arg, typename Handler>
  //   static Gen::type_normal<Arg, Handler>::type
  //   Gen::generate_normal(const Gen&, const transport&);
  //
  //   template <typename Arg, typename Handler>
  //   struct Gen::type_duplicate_removal {
  //     typename type;
  //     requires MessageType<type>;
  //   };
  //   template <typename Arg, typename Handler>
  //   static Gen::type_duplicate_removal<Arg, Handler>::type
  //   Gen::generate_duplicate_removal(const Gen&, const transport&);
  //
  //   template <typename Arg, typename Handler, typename Less>
  //   struct Gen::type_suboptimal_value_removal {
  //     typename type;
  //     requires MessageType<type>;
  //   };
  //   template <typename Arg, typename Handler, typename Less>
  //   static Gen::type_suboptimal_value_removal<Arg, Handler, Less>::type
  //   Gen::generate_suboptimal_value_removal(const Gen&, const transport&, const Less&);
  //
  //   template <typename Arg, typename Handler, typename Combine>
  //   struct Gen::type_combination {
  //     typename type;
  //     requires MessageType<type>;
  //   };
  //   template <typename Arg, typename Handler, typename Combine>
  //   static Gen::type_combination<Arg, Handler, Combine>::type
  //   Gen::generate_combination(const Gen&, const transport&, const Combine&);
  // }
  
  BOOST_PARAMETER_NAME(generator)
  BOOST_PARAMETER_NAME(transport)
  BOOST_PARAMETER_NAME(arg_type)
  BOOST_PARAMETER_NAME(handler_type)
  BOOST_PARAMETER_NAME(owner)
  BOOST_PARAMETER_NAME(sent_part)
  BOOST_PARAMETER_NAME(reduction)

  enum no_reduction_t {no_reduction};

  struct make_default_policy {};
  template <typename F, typename Arg> struct duplicate_policy_from_projection;
  template <typename Project, typename Policy, typename Arg>
  struct get_policy_type {typedef Policy type; static type get_policy(const Project&, const Policy& p) {return p;}};
  template <typename Project, typename Arg>
  struct get_policy_type<Project, make_default_policy, Arg> {typedef duplicate_policy_from_projection<Project, Arg> type; static type get_policy(const Project& pr, const make_default_policy&) {return type(pr);}};

  template <typename Project, typename Policy = make_default_policy>
  struct duplicate_removal_t {
    Project project;
    Policy policy;
    explicit duplicate_removal_t(const Project& p): project(p), policy() {}
    duplicate_removal_t(const Project& p, const Policy& po): project(p), policy(po) {}
  };

  template <typename Less>
  struct suboptimal_value_removal_t {Less less; suboptimal_value_removal_t(const Less& less): less(less) {}};

  template <typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
  struct combination_t {
    Combine combine;
    OptIdentity opt_identity;
    GetKey get_key;
    GetValue get_value;
    MakeKeyval make_keyval;
    combination_t(const Combine& combine, const OptIdentity& opt_identity, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval)
      : combine(combine), opt_identity(opt_identity), get_key(get_key), get_value(get_value), make_keyval(make_keyval) {}
  };
  template <typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
  struct idempotent_combination_t {
    Combine combine;
    OptIdentity opt_identity;
    GetKey get_key;
    GetValue get_value;
    MakeKeyval make_keyval;
    idempotent_combination_t(const Combine& combine, const OptIdentity& opt_identity, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval)
      : combine(combine), opt_identity(opt_identity), get_key(get_key), get_value(get_value), make_keyval(make_keyval) {}
  };

  template <typename Project> duplicate_removal_t<Project> duplicate_removal(const Project& p) {return duplicate_removal_t<Project>(p);}
  template <typename Project, typename Policy> duplicate_removal_t<Project, Policy> duplicate_removal(const Project& p, const Policy& po) {return duplicate_removal_t<Project, Policy>(p, po);}

  template <typename Combine>
  combination_t<Combine, no_identity_t, pair_first, pair_second, make_pair_t>
  combination(const Combine& c) {
    return combination_t<Combine, no_identity_t, pair_first, pair_second, make_pair_t>(c, no_identity_t(), pair_first(), pair_second(), make_pair_t());
  }

  template <typename Combine, typename Identity>
  combination_t<Combine, Identity, pair_first, pair_second, make_pair_t>
  combination(const Combine& c, const Identity& i) {
    return combination_t<Combine, Identity, pair_first, pair_second, make_pair_t>(c, i, pair_first(), pair_second(), make_pair_t());
  }

  template <typename Combine, typename Identity, typename GetKey, typename GetValue, typename MakeKeyval>
  combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval>
  combination(const Combine& c, const Identity& i, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval) {
    return combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval>(c, i, get_key, get_value, make_keyval);
  }

  template <typename Combine>
  idempotent_combination_t<Combine, no_identity_t, pair_first, pair_second, make_pair_t>
  idempotent_combination(const Combine& c) {
    return idempotent_combination_t<Combine, no_identity_t, pair_first, pair_second, make_pair_t>(c, no_identity_t(), pair_first(), pair_second(), make_pair_t());
  }

  template <typename Combine, typename Identity>
  idempotent_combination_t<Combine, Identity, pair_first, pair_second, make_pair_t>
  idempotent_combination(const Combine& c, const Identity& i) {
    return idempotent_combination_t<Combine, Identity, pair_first, pair_second, make_pair_t>(c, i, pair_first(), pair_second(), make_pair_t());
  }

  template <typename Combine, typename Identity, typename GetKey, typename GetValue, typename MakeKeyval>
  idempotent_combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval>
  idempotent_combination(const Combine& c, const Identity& i, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval) {
    return idempotent_combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval>(c, i, get_key, get_value, make_keyval);
  }

  template <typename> struct is_valid_reduction: boost::mpl::false_ {};
  template <> struct is_valid_reduction<no_reduction_t>: boost::mpl::true_ {};
  template <typename Project, typename Policy> struct is_valid_reduction<duplicate_removal_t<Project, Policy> >: boost::mpl::true_ {};
  template <typename Combine, typename Identity, typename GetKey, typename GetValue, typename MakeKeyval>
  struct is_valid_reduction<combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval> >: boost::mpl::true_ {};
  template <typename Combine, typename Identity, typename GetKey, typename GetValue, typename MakeKeyval>
  struct is_valid_reduction<idempotent_combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval> >: boost::mpl::true_ {};

  // Simpler interface
  template <typename ArgType, typename HandlerType, typename Gen, typename Owner, typename Reduction = no_reduction_t>
  struct message_type_generator_result
    : Gen::template call_result<ArgType, HandlerType, Owner, Reduction> {};

  struct noncoalesced_generator {
    template <typename Arg, typename Handler, typename Owner, typename Reduction = no_reduction_t>
    struct call_result {
      struct type {
        Owner owner;
        message_type<Arg> mt;
        boost::shared_ptr<Handler> h;

        type(const noncoalesced_generator&, transport trans, const Owner& owner, const Reduction& = no_reduction)
          : owner(owner), mt(trans.create_message_type<Arg>()), h()
        {
          mt.set_max_count(1);
          mt.set_handler(boost::function<void (transport::rank_type, const Arg*, size_t)>(boost::ref(*this)));
        }
        struct eat_sp {void operator()(const boost::shared_ptr<Arg>&) const {}};
        void send(const Arg& a) {
          transport::rank_type dest = get(owner, a);
          mt.message_being_built(dest);
          boost::shared_ptr<Arg> buf(new Arg(a));
          mt.send(buf.get(), 1, dest, boost::bind<void>(eat_sp(), buf)); // Keep ownership of buf
        }
        void send_with_tid(const Arg& a, int /*tid*/) {
          this->send(a);
        }
        void set_handler(const Handler& h_) {h.reset(new Handler(h_));}
#ifndef BOOST_NO_RVALUE_REFERENCES
        void set_handler(Handler&& h_) {h.reset(new Handler(std::move(h_)));}
#endif
        const Handler& get_handler() const {assert (h); return *h;}

        typedef void result_type;
        void operator()(transport::rank_type, const Arg* a, size_t count) const {
          assert (count == 1); (void)count;
          assert (h);
          (*h)(*a);
        }
      };
    };
  };

  

  template <typename CoalescingGen>
  struct simple_generator {
    CoalescingGen coalescing_gen;

    explicit simple_generator(const CoalescingGen& coalescing_gen)
      : coalescing_gen(coalescing_gen)
    {}

    template <typename Arg, typename Handler, typename Owner, typename Reduction = no_reduction_t, typename BufferSorter = DummyBufferSorter<Arg> >
    struct call_result {
      struct type: object_based_addressing<Arg, Handler, Owner, CoalescingGen, BufferSorter> {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
	detail::hit_rate_counter counters;
#endif
        type(const simple_generator& gen, transport trans, const Owner& owner, const Reduction& = no_reduction, const BufferSorter& bufsrter = DummyBufferSorter<Arg>())
          : 
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
	  counters(trans.get_nthreads()),
#endif
	  object_based_addressing<Arg, Handler, Owner, CoalescingGen, BufferSorter>(gen.coalescing_gen, trans, owner, bufsrter) {} // need buffer sorter here
      };
    };
  };

  template <typename CoalescingGen, typename Routing>
  struct routing_generator {
    CoalescingGen cg;
    Routing routing;

    explicit routing_generator(const CoalescingGen& c, const Routing& routing): cg(c), routing(routing) {}

    template <typename Arg, typename Handler, typename Owner, typename Reduction = no_reduction_t, typename BufferSorter = DummyBufferSorter<Arg> >
    struct call_result {
      struct type: object_based_addressing_dest_hbr<Arg, Handler, Owner, Routing, CoalescingGen, BufferSorter> {
	#ifdef AMPLUSPLUS_PRINT_HIT_RATES
	detail::hit_rate_counter counters;
#endif
        type(const routing_generator& gen, amplusplus::transport trans, const Owner& owner, const Reduction& = no_reduction, const BufferSorter& bufsrter = DummyBufferSorter<Arg>())
          : 
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
	  counters(trans.get_nthreads()),
#endif
	  object_based_addressing_dest_hbr<Arg, Handler, Owner, Routing, CoalescingGen, BufferSorter>(std::make_pair(gen.cg, gen.routing), trans, owner, bufsrter) {}
      };
    };
  };

  template <typename F, typename Arg>
  struct duplicate_policy_from_projection {
    typedef Arg value_type;
    typedef typename boost::result_of<F(Arg)>::type stored_type;

    private:
    F f;

    public:
    explicit duplicate_policy_from_projection(const F& f = F()): f(f) {}

    stored_type stored_value(const Arg& d) const {
      return f(d);
    }
    bool match(const Arg& d, stored_type stored) const {
      return f(d) == stored;
    }
    size_t hash(stored_type v, size_t size) const {
      size_t key = boost::hash<stored_type>()(v);
      key ^= (key / size);
      key ^= (key / size / size);
      return key %= size;
    }
    stored_type dummy_value(size_t h, size_t size) const {
      for(size_t i = 0; ; ++i) {
	stored_type x(h ^ i);
	if(hash(x, size) != h)
	  return x;
      }
      return stored_type(); // We should never be here
    }
  };

  namespace detail {
    template <typename Arg, typename Handler, typename Owner, typename Routing, typename CoalescingGen, typename BufferSorter>
    struct make_oba {
      struct type: object_based_addressing_dest_hbr<Arg, Handler, Owner, Routing, CoalescingGen, BufferSorter> {
	template <typename Gen>
	type(const Gen& gen, transport trans, const Owner& owner, const Routing& r, const BufferSorter& bufstr):
	  object_based_addressing_dest_hbr<Arg, Handler, Owner, Routing, CoalescingGen, BufferSorter>
          (std::make_pair(gen, r), trans, owner, bufstr) {}
      };
    };
    
    template <typename Arg, typename Handler, typename Owner, typename CoalescingGen, typename BufferSorter>
    struct make_oba<Arg, Handler, Owner, no_routing, CoalescingGen, BufferSorter> {
      struct type: object_based_addressing<Arg, Handler, Owner, CoalescingGen, BufferSorter> {
	template <typename Gen>
	type(const Gen& gen, transport trans, const Owner& owner, const no_routing&, const BufferSorter& bufstr):
	  object_based_addressing<Arg, Handler, Owner, CoalescingGen, BufferSorter>(gen, trans, owner, bufstr) {}
      };
    };
  }

  namespace detail {
    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Reduction, typename BufferSorter = DummyBufferSorter<Arg> >
    struct cache_generator_call_result;
  }
    
  template <typename CoalescingGen, typename Routing>
  struct cache_generator {
    CoalescingGen cg;
    Routing routing;
    unsigned int lg_size;

    explicit cache_generator(const CoalescingGen& c, unsigned int lg_size, const Routing& routing = no_routing())
      : cg(c), routing(routing), lg_size(lg_size) {}

    template <typename Arg, typename Handler, typename Owner, typename Reduction, typename BufferSorter = DummyBufferSorter<Arg> >
    struct call_result: detail::cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, Reduction, BufferSorter> {};
  };

  namespace detail {
    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Reduction, typename BufferSorter>
    struct cache_generator_call_result {
      typedef cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, CoalescingGen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const Reduction& r, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(gen.cg, trans, owner, gen.routing, bufstr) {}
      };
    };

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Project, typename Policy, typename BufferSorter>
    struct cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, duplicate_removal_t<Project, Policy>, BufferSorter> {
      typedef amplusplus::simple_cache_remove_duplicates_gen<CoalescingGen, typename get_policy_type<Project, Policy, Arg>::type> message_type_gen;
      typedef cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, message_type_gen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const duplicate_removal_t<Project, Policy>& r, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(message_type_gen(gen.cg, get_policy_type<Project, Policy, Arg>::get_policy(r.project, r.policy), gen.lg_size), trans, owner, gen.routing, bufstr) {}
      };
    };

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Combine, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval, typename BufferSorter>
    struct cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>, BufferSorter> {
      typedef amplusplus::locking_cache_binop_reduction_gen<CoalescingGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval> message_type_gen;
      typedef cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, message_type_gen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>& comb, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(message_type_gen(gen.cg, comb.combine, comb.opt_identity, comb.get_key, comb.get_value, comb.make_keyval, gen.lg_size), trans, owner, gen.routing, bufstr) {}
      };
    };

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Combine, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval, typename BufferSorter>
    struct cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, idempotent_combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>, BufferSorter>
      : cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>, BufferSorter> {};

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Reduction, typename BufferSorter = DummyBufferSorter<Arg> >
    struct per_thread_cache_generator_call_result;
  } // namespace detail

  template <typename CoalescingGen, typename Routing>
  struct per_thread_cache_generator {
    CoalescingGen cg;
    Routing routing;
    unsigned int lg_size;
    
    explicit per_thread_cache_generator(const CoalescingGen& c, unsigned int lg_size, const Routing& routing = no_routing())
      : cg(c), routing(routing), lg_size(lg_size) {}
    
    template <typename Arg, typename Handler, typename Owner, typename Reduction, typename BufferSorter = DummyBufferSorter<Arg> >
    struct call_result: detail::per_thread_cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, Reduction, BufferSorter> {};
  };
  
  namespace detail {
    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Reduction, typename BufferSorter>
    struct per_thread_cache_generator_call_result {
      typedef per_thread_cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, CoalescingGen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const Reduction& r, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(gen.cg, trans, owner, gen.routing, bufstr) {}
      };
    };

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Project, typename Policy, typename BufferSorter>
    struct per_thread_cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, duplicate_removal_t<Project, Policy>, BufferSorter> {
      typedef amplusplus::per_thread_cache_remove_duplicates_gen<CoalescingGen, typename get_policy_type<Project, Policy, Arg>::type> message_type_gen;
      typedef per_thread_cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, message_type_gen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const duplicate_removal_t<Project, Policy>& r, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(message_type_gen(gen.cg, get_policy_type<Project, Policy, Arg>::get_policy(r.project, r.policy), gen.lg_size), trans, owner, gen.routing, bufstr) {}
      };
    };

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Combine, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval, typename BufferSorter>
    struct per_thread_cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>, BufferSorter> {
      typedef amplusplus::per_thread_cache_binop_reduction_gen<CoalescingGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval> message_type_gen;
      typedef per_thread_cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, message_type_gen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>& comb, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(message_type_gen(gen.cg, comb.combine, comb.opt_identity, comb.get_key, comb.get_value, comb.make_keyval, gen.lg_size), owner, gen.routing, bufstr) {}
      };
    };

    template <typename CoalescingGen, typename Routing, typename Arg, typename Handler, typename Owner, typename Combine, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval, typename BufferSorter>
    struct per_thread_cache_generator_call_result<CoalescingGen, Routing, Arg, Handler, Owner, idempotent_combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>, BufferSorter > {
      typedef amplusplus::per_thread_wt_cache_binop_reduction_gen<CoalescingGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval> message_type_gen;
      typedef per_thread_cache_generator<CoalescingGen, Routing> cgen;
      typedef detail::make_oba<Arg, Handler, Owner, Routing, message_type_gen, BufferSorter> type_base;
      struct type: type_base::type {
	type(const cgen& gen, transport trans, const Owner& owner, const idempotent_combination_t<Combine, OptIdentity, GetKey, GetValue, MakeKeyval>& comb, const BufferSorter& bufstr = DummyBufferSorter<Arg>())
	  : type_base::type(message_type_gen(gen.cg, comb.combine, comb.opt_identity, comb.get_key, comb.get_value, comb.make_keyval, gen.lg_size), trans, owner, gen.routing, bufstr) {}
      };
    };
  } // namespace detail
}

#endif // AMPLUSPLUS_MESSAGE_TYPE_GENERATORS_HPP
