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

#ifndef AMPLUSPLUS_REDUCTIONS_HPP
#define AMPLUSPLUS_REDUCTIONS_HPP

#include <am++/traits.hpp>
#include <boost/assert.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/unordered_set.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <boost/serialization/static_warning.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/int.hpp>
#include <vector>
#include <iostream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <am++/detail/thread_support.hpp>
#include <am++/detail/signal.hpp>
#include <am++/detail/factory_wrapper.hpp>
#include <am++/dummy_buffer_sorter.hpp>

// #define AMPLUSPLUS_PRINT_HIT_RATES

namespace amplusplus {

#if 0
concept DuplicateRemovalPolicy<typename P> {
  typename P::value_type; // Message type
  typename P::stored_type; // Must be able to be used with amplusplus::detail::atomic
  stored_type P::stored_value(value_type) const; // Get value to put into table
  bool P::match(value_type a, stored_type b) const; // a == b, match(x, stored_value(x)) must be true
  size_t P::hash(stored_type) const; // Generate a hash
  stored_type P::dummy_value(size_t h) const; // Generate a value that either will not be in the input or does not have hash code h
}

concept SuboptimalValueRemovalPolicy<typename P> {
  typename P::value_type; // Message type
  void P::combine(value_type a, value_type b, value_type& combined_for_new_cache, bool& send_a, bool& send_b) const; // a is cached value, b is new value, replace a by combined_for_new_cache, if send_a then send a, if send_b then send b
  size_t P::hash(value_type) const; // Generate a hash
  value_type P::dummy_value(size_t h) const; // Generate a value that either will not be in the input or does not have hash code h
  bool P::is_valid(size_t h, value_type v) const; // Test whether v is a valid value for hash h
}

concept BinaryOpPolicy<typename P> {
  typename P::pair_type;
  typename P::key_type;
  typename P::value_type;
  key_type& P::get_key(pair_type&) const;
  const key_type& P::get_key(const pair_type&) const;
  value_type& P::get_value(pair_type&) const;
  const value_type& P::get_value(const pair_type&) const;
  pair_type P::make_pair(key_type, value_type) const;
  size_t P::hash_key(key_type) const;
  key_type P::dummy_key(size_t h) const;
  bool P::is_valid_key(size_t h, key_type) const;
  value_type P::op(value_type, value_type) const;
  bool P::is_identity(value_type) const;
}
#endif

  namespace detail {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
    struct hit_rate_counter {
      std::vector<size_t> hits, tests;
      hit_rate_counter(unsigned int no_of_threads): hits(no_of_threads, 0), tests(no_of_threads, 0) {}
      void clear() {
	for(unsigned int i = 0; i < hits.size(); ++i) {
	  hits[i] = 0;
	  tests[i] = 0;
	}
      }
      void test(unsigned int tid) { ++tests[tid]; }
      void hit(unsigned int tid) { ++hits[tid]; }
#if 0
      ~hit_rate_counter() {
        std::ostringstream os;
        os << "Tests = " << tests.load() << ", hits = " << hits.load() << ", hit rate = " << (double(hits.load()) / double(tests.load()) * 100.) << "%" << std::endl;
        std::cerr << os.str() << std::flush;
      }
#endif
    };
#else
    struct hit_rate_counter {
      hit_rate_counter(unsigned int no_of_threads) {}
      void clear() {}
      void test(unsigned int tid) {}
      void hit(unsigned int tid) {}
    };
#endif

    template <size_t N>
    struct is_power_of_2: boost::mpl::and_<boost::mpl::bool_<N % 2 == 0>, is_power_of_2<N / 2> > {};
    template <>
    struct is_power_of_2<1>: boost::mpl::true_ {};
  }

template <typename CoalescingLayer, typename Policy>
class simple_cache_remove_duplicates;

template <typename CoalescingLayerGen, typename Policy>
struct simple_cache_remove_duplicates_gen {
  CoalescingLayerGen clg;
  Policy policy;
  unsigned int lg_size;
  simple_cache_remove_duplicates_gen(const CoalescingLayerGen& clg, const Policy& policy, unsigned int lg_size): clg(clg), policy(policy), lg_size(lg_size) {}

  template <typename Arg, typename Handler>
  struct inner {
    typedef simple_cache_remove_duplicates<typename CoalescingLayerGen::template inner<Arg, Handler>::type, Policy> type;
  };
};

template <typename CoalescingLayer, typename Policy>
class simple_cache_remove_duplicates: public CoalescingLayer {
  public:
  typedef typename message_type_traits<CoalescingLayer>::arg_type arg_type;
  typedef typename message_type_traits<CoalescingLayer>::handler_type handler_type;

  private:
  typedef transport::rank_type rank_type;
  typedef typename Policy::stored_type stored_type;

  Policy policy;
  unsigned int lg_size;
public:
  detail::hit_rate_counter counters;
private:
  valid_rank_set possible_dests;
  transport::rank_type num_possible_dests;
  boost::scoped_array<amplusplus::detail::atomic<stored_type> > cache; // possible_dests->count() chunks of this->get_size() elements each
  boost::scoped_array<amplusplus::detail::atomic<stored_type>*> start_ptrs; // Beginning of cache for rank or NULL

  public:
  template <typename CoalescingLayerGen>
  explicit simple_cache_remove_duplicates(
             const simple_cache_remove_duplicates_gen<CoalescingLayerGen, Policy>& gen,
             const transport& trans,
             const valid_rank_set& possible_dests,
             const valid_rank_set& /*possible_sources*/)
      : CoalescingLayer(gen.clg, trans),
        policy(gen.policy),
        lg_size(gen.lg_size),
        counters(trans.get_nthreads()),
        possible_dests(possible_dests),
        num_possible_dests(possible_dests->count()),
        cache(new amplusplus::detail::atomic<stored_type>[possible_dests->count() * (size_t(1) << this->lg_size)]),
        start_ptrs(new amplusplus::detail::atomic<stored_type>*[this->get_transport().size()])
    {
      std::fill(&start_ptrs[0], &start_ptrs[this->get_transport().size()], (amplusplus::detail::atomic<stored_type>*)0);
      for (transport::rank_type i = 0; i < num_possible_dests; ++i) {
        start_ptrs[possible_dests->rank_from_index(i)] = &cache[i * (size_t(1) << this->lg_size)];
      }
      clear();
    }

  size_t get_size() const {return size_t(1) << lg_size;}

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  private: simple_cache_remove_duplicates(const simple_cache_remove_duplicates&); public:
#else
  simple_cache_remove_duplicates(const simple_cache_remove_duplicates&) = delete;
#endif

  void clear() {
    counters.clear();
    for (size_t i = 0; i < num_possible_dests * this->get_size(); ++i) {
      // Ensure this entry can never be found to avoid false positives
      cache[i].store(policy.dummy_value(i % this->get_size(), this->get_size()));
    }
  }

  void send(const arg_type& arg, rank_type dest) {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.test(amplusplus::detail::get_thread_id());
#endif
    BOOST_ASSERT (this->get_transport().is_valid_rank(dest));
    stored_type new_val = policy.stored_value(arg);
    size_t key = policy.hash(new_val, this->get_size());
    amplusplus::detail::atomic<stored_type>* start_ptr = start_ptrs[dest];
    BOOST_ASSERT (start_ptr);
    stored_type old_val = start_ptr[key].exchange(new_val);
    if (!policy.match(arg, old_val)) {
      CoalescingLayer::send(arg, dest);
    } else {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.hit(amplusplus::detail::get_thread_id());
#endif
    }
  }

  void send_with_tid(const arg_type& arg, rank_type dest, int tid) {
    counters.test(tid);
    BOOST_ASSERT (this->get_transport().is_valid_rank(dest));
    stored_type new_val = policy.stored_value(arg);
    size_t key = policy.hash(new_val, this->get_size());
    amplusplus::detail::atomic<stored_type>* start_ptr = start_ptrs[dest];
    BOOST_ASSERT (start_ptr);
    stored_type old_val = start_ptr[key].exchange(new_val);
    if (!policy.match(arg, old_val)) {
      CoalescingLayer::send_with_tid(arg, dest, tid);
    } else {
      counters.hit(tid);
    }
  }

  bool flush(transport::rank_type r) {clear(); return false;}

  const handler_type& get_handler() {return CoalescingLayer::get_handler();}
  void set_handler(const handler_type& h) {CoalescingLayer::set_handler(h);}
  transport get_transport() const {return CoalescingLayer::get_transport();}
};

template <typename CoalescingLayer, typename Policy, typename BufferSorter>
class per_thread_cache_remove_duplicates;

template <typename CoalescingLayerGen, typename Policy>
struct per_thread_cache_remove_duplicates_gen {
  CoalescingLayerGen clg;
  Policy policy;
  unsigned int lg_size;
  per_thread_cache_remove_duplicates_gen(const CoalescingLayerGen& clg, const Policy& policy, unsigned int lg_size): clg(clg), policy(policy), lg_size(lg_size) {}

  template <typename Arg, typename Handler, typename BufferSorter = DummyBufferSorter<Arg> >
  struct inner {
    typedef per_thread_cache_remove_duplicates<typename CoalescingLayerGen::template inner<Arg, Handler, BufferSorter>::type, Policy, BufferSorter> type;
  };
};

template <typename CoalescingLayer, typename Policy, typename BufferSorter>
class per_thread_cache_remove_duplicates: public CoalescingLayer {
  public:
  typedef typename message_type_traits<CoalescingLayer>::arg_type arg_type;
  typedef typename message_type_traits<CoalescingLayer>::handler_type handler_type;

  private:
  typedef transport::rank_type rank_type;
  typedef typename Policy::stored_type stored_type;

  Policy policy;
  unsigned int lg_size;
public:
  detail::hit_rate_counter counters;
private:
  valid_rank_set possible_dests;
  transport::rank_type num_possible_dests;
  transport::rank_type num_ranks;
  int num_threads;
  boost::scoped_array<stored_type> cache; // possible_dests->count() chunks of this->get_size() elements each
  boost::scoped_array<stored_type*> start_ptrs; // Beginning of cache for rank or NULL

  public:
  template <typename CoalescingLayerGen>
  explicit per_thread_cache_remove_duplicates(
             const per_thread_cache_remove_duplicates_gen<CoalescingLayerGen, Policy>& gen,
             const transport& trans,
             const valid_rank_set& possible_dests,
             const valid_rank_set& possible_sources,
             const BufferSorter bufsrter
             )
      : CoalescingLayer(gen.clg, trans, possible_dests, possible_sources, bufsrter),
        policy(gen.policy),
        lg_size(gen.lg_size),
        counters(trans.get_nthreads()),
        possible_dests(possible_dests),
        num_possible_dests(possible_dests->count()),
        num_ranks(this->get_transport().size()),
        num_threads(this->get_transport().get_nthreads()),
        cache(new stored_type[this->num_threads * this->num_possible_dests * (size_t(1) << this->lg_size)]),
        start_ptrs(new stored_type*[this->num_threads * this->num_ranks])
    {
      std::fill(&start_ptrs[0],
                &start_ptrs[this->num_threads * this->num_ranks],
                (stored_type*)0);
      for (transport::rank_type i = 0; i < num_possible_dests * num_threads; ++i) {
        BOOST_ASSERT (possible_dests->rank_from_index(i % num_possible_dests) < num_ranks);
        start_ptrs[(i / num_possible_dests) * num_ranks + possible_dests->rank_from_index(i % num_possible_dests)]
          = &cache[i * this->get_size()];
      }
      clear();
    }

  size_t get_size() const {return size_t(1) << lg_size;}

  void clear() {
    counters.clear();
    for (size_t i = 0; i < num_threads * num_possible_dests * this->get_size(); ++i) {
      // Ensure this entry can never be found to avoid false positives
      cache[i] = policy.dummy_value(i % this->get_size(), this->get_size());
      BOOST_ASSERT (policy.hash(cache[i], this->get_size()) != i % this->get_size());
    }
  }

  void send_with_tid(const arg_type& arg, rank_type dest, const int thread_id) {
    BOOST_ASSERT (thread_id >= 0 && thread_id < (int)num_threads);
    counters.test(thread_id);
    BOOST_ASSERT (this->get_transport().is_valid_rank(dest));
    BOOST_ASSERT (dest < num_ranks);
    stored_type new_val = policy.stored_value(arg);
    size_t h = policy.hash(new_val, this->get_size());
    BOOST_ASSERT (h < this->get_size());
    stored_type* start_ptr = start_ptrs[(size_t)thread_id * num_ranks + dest];
    BOOST_ASSERT (start_ptr);
    if (!policy.match(arg, start_ptr[h])) {
      start_ptr[h] = new_val;
      CoalescingLayer::send_with_tid(arg, dest, thread_id);
    } else {
      counters.hit(thread_id);
    }
  }

  void send(const arg_type& arg, rank_type dest) {
    this->send_with_tid(arg, dest, amplusplus::detail::get_thread_id());
  }

  bool flush(transport::rank_type r) {
    const int thread_id = amplusplus::detail::get_thread_id();
    BOOST_ASSERT (thread_id >= 0 && thread_id < (int)num_threads);
    BOOST_ASSERT (r < num_ranks);
    stored_type* start_ptr = start_ptrs[(size_t)thread_id * num_ranks + r];
    if (start_ptr != 0) { // Skip if not a valid destination
      for (size_t i = 0; i < this->get_size(); ++i) {
        start_ptr[i] = policy.dummy_value(i);
        BOOST_ASSERT (policy.hash(start_ptr[i], this->get_size()) != i % this->get_size());
      }
    }
    return false;
  }

  const handler_type& get_handler() {return CoalescingLayer::get_handler();}
  void set_handler(const handler_type& h) {CoalescingLayer::set_handler(h);}
#ifndef BOOST_NO_RVALUE_REFERENCES
  void set_handler(handler_type&& h) {CoalescingLayer::set_handler(std::move(h));}
#endif
  transport get_transport() const {return CoalescingLayer::get_transport();}
};

template <typename CoalescingLayerGen, typename Policy, typename Arg, typename Handler>
class unordered_set_remove_duplicates;

template <typename CoalescingLayerGen, typename Policy>
struct unordered_set_remove_duplicates_gen {
  CoalescingLayerGen clg;
  Policy policy;
  unordered_set_remove_duplicates_gen(const CoalescingLayerGen& clg, const Policy& policy): clg(clg), policy(policy) {}

  template <typename Arg, typename Handler>
  struct inner {
    typedef unordered_set_remove_duplicates<CoalescingLayerGen, Policy, Arg, Handler> type;
  };
};

template <typename CoalescingLayerGen, typename Policy, typename Arg, typename Handler>
class unordered_set_remove_duplicates {
  struct policy_hasher {
    const Policy& p;
    policy_hasher(const Policy& p): p(p) {}
    size_t operator()(const typename Policy::stored_type& v) const {return p.hash(v);}
  };

  public:
  typedef Arg arg_type;
  typedef Handler handler_type;

  private:
  typedef transport::rank_type rank_type;
  typedef typename Policy::stored_type stored_type;

  typename CoalescingLayerGen::template inner<arg_type, handler_type>::type cl;
  Policy policy;
public:
  detail::hit_rate_counter counters;
private:
  typedef boost::unordered_set<stored_type, policy_hasher> uo_set_type;
  std::vector<uo_set_type> cache; // One set for each destination, vector rather than scoped_array because our hasher is not default constructible
  boost::scoped_array<boost::mutex> lock; // One for each destination

  public:
  explicit unordered_set_remove_duplicates(
             const unordered_set_remove_duplicates_gen<CoalescingLayerGen, Policy>& gen,
             const transport& trans,
             const valid_rank_set& possible_dests,
             const valid_rank_set& possible_sources)
      : cl(gen.clg, trans, possible_dests, possible_sources),
        policy(gen.policy),
        counters(trans.get_nthreads()),
        cache(cl.get_transport().size(), uo_set_type(100, policy_hasher(this->policy))),
        lock(AMPLUSPLUS_MULTITHREAD(new boost::mutex[cl.get_transport().size()]))
    {clear();}

  void clear() {
    counters.clear();
    rank_type nranks = cl.get_transport().size();
    for (rank_type i = 0; i < nranks; ++i) {
      AMPLUSPLUS_MULTITHREAD(boost::lock_guard<boost::mutex> l(lock[i]);)
      cache[i].clear();
    }
  }

  void send(const arg_type& arg, rank_type dest) {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
    counters.test(amplusplus::detail::get_thread_id());
#endif
    BOOST_ASSERT (cl.get_transport().is_valid_rank(dest));
    stored_type new_val = policy.stored_value(arg);
    AMPLUSPLUS_MULTITHREAD(boost::lock_guard<boost::mutex> l(lock[dest]);)
    if (cache[dest].insert(new_val).second) { // If not present before, insert and return true
      cl.send(arg, dest);
      if (cache[dest].size() > 10000) cache[dest].clear();
    } else {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.hit(amplusplus::detail::get_thread_id());
#endif
    }
  }

  void send_with_tid(const arg_type& arg, rank_type dest, int tid) {
    counters.test(tid);
    BOOST_ASSERT (cl.get_transport().is_valid_rank(dest));
    stored_type new_val = policy.stored_value(arg);
    AMPLUSPLUS_MULTITHREAD(boost::lock_guard<boost::mutex> l(lock[dest]);)
    if (cache[dest].insert(new_val).second) { // If not present before, insert and return true
      cl.send_with_tid(arg, dest, tid);
      if (cache[dest].size() > 10000) cache[dest].clear();
    } else {
      counters.hit(tid);
    }
  }

};

struct no_identity_t {};

template <typename OptIdentity, typename GetValue, typename Arg>
bool is_identity(const OptIdentity& opt_identity, const GetValue& get_value, const Arg& arg) {
  return get_value(arg) == opt_identity;
}

template <typename GetValue, typename Arg>
bool is_identity(const no_identity_t&, const GetValue&, const Arg&) {return false;}

struct pair_first {
  template <typename P>
  typename P::first_type operator()(const P& p) const {return p.first;}

  template <typename F>
  struct result {
    typedef typename boost::function_traits<F>::arg1_type::first_type type;
  };
};

struct pair_second {
  template <typename P>
  typename P::second_type operator()(const P& p) const {return p.second;}

  template <typename F>
  struct result {
    typedef typename boost::function_traits<F>::arg1_type::second_type type;
  };
};

struct make_pair_t {
  template <typename A, typename B>
  std::pair<A, B> operator()(const A& a, const B& b) const {return std::make_pair(a, b);}

  template <typename F>
  struct result {
    typedef std::pair<typename boost::function_traits<F>::arg1_type, typename boost::function_traits<F>::arg2_type> type;
  };
};

template <typename T, typename = void> struct dummy_value;

template <typename T>
struct dummy_value<T, typename boost::enable_if<boost::is_integral<T> >::type> {
  T operator()(size_t h) const {return T(h ^ 1);}
};

template <typename CoalescingLayer, typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
class simple_cache_binop_reduction;

template <typename CoalescingLayer, typename Combine, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval>
class simple_cache_binop_reduction {
  public:
  typedef typename message_type_traits<CoalescingLayer>::arg_type arg_type;
  typedef typename message_type_traits<CoalescingLayer>::handler_type handler_type;

  private:
  typedef typename boost::result_of<GetKey(arg_type)>::type key_type;
  typedef typename boost::result_of<GetValue(arg_type)>::type value_type;

  typedef transport::rank_type rank_type;

  CoalescingLayer& cl;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  boost::hash<key_type> hash_key;
  unsigned int lg_size;
  boost::scoped_array<amplusplus::detail::atomic<arg_type> > values;

  public:
  detail::hit_rate_counter counters;

  explicit simple_cache_binop_reduction(CoalescingLayer& cl, unsigned int lg_size, const Combine& combine = Combine(), const OptIdentity& opt_identity = OptIdentity(), const GetKey& get_key = GetKey(), const GetValue& get_value = GetValue(), const MakeKeyval& make_keyval = MakeKeyval())
    : cl(cl),
      combine(combine),
      opt_identity(opt_identity),
      get_key(get_key),
      get_value(get_value),
      make_keyval(make_keyval),
      hash_key(),
      lg_size(lg_size),
      values(new amplusplus::detail::atomic<arg_type>[(size_t(1) << this->lg_size) * cl.get_transport().size()]),
      counters(cl.get_transport().get_nthreads())
  {
    this->clear();
  }

  size_t get_size() const {return size_t(1) << lg_size;}

  void clear() { // Clear without sending data
    counters.clear();
    size_t nranks = cl.get_transport().size();
    for (size_t i = 0; i < this->get_size() * nranks; ++i) {
      values[i].store(make_keyval(get_dummy_value(i % this->get_size()), value_type()));
    }
  }

  bool flush(rank_type dest) { // Clear while sending any valid data
    bool did_anything = false;
    for (size_t h = 0; h < this->get_size(); ++h) {
      size_t index = dest * this->get_size() + h;
      arg_type old_kv = values[index].exchange(make_keyval(get_dummy_value(h), value_type()));
      if (hash_key(get_key(old_kv)) == h && !is_identity(opt_identity, get_value, old_kv)) {
        cl.send(old_kv, dest);
        did_anything = true;
      }
    }
    return did_anything;
  }

  void send(const arg_type& a, rank_type dest) {
    BOOST_STATIC_ASSERT (boost::is_void<arg_type>::value && !"Fix me up for new flush framework");
    key_type a_key = get_key(a);
    value_type a_value = get_value(a);
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.test(amplusplus::detail::get_thread_id());
#endif
    size_t h = hash_key(a_key);
    h ^= (h / this->get_size());
    h ^= (h / this->get_size() / this->get_size());
    h %= this->get_size();
    size_t index = dest * this->get_size() + h;
    arg_type old_kv = values[index].load();
    while (true) {
      key_type old_key = get_key(old_kv);
      value_type old_value = get_value(old_kv);
      if (a_key == old_key) {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.hit(amplusplus::detail::get_thread_id());
#endif
        value_type new_value = combine(old_value, a_value);
        arg_type new_kv = make_keyval(a_key, new_value);
        if (new_value == old_value) break;
        bool changed = values[index].compare_exchange_strong(old_kv, new_kv);
        // When changed is false, compare_exchange_strong updated old_kv, so go
        // to top of loop with that value
        if (changed) break; else continue;
      } else {
        bool changed = values[index].compare_exchange_strong(old_kv, a);
        // When changed is false, compare_exchange_strong updated old_kv, so go
        // to top of loop with that value
        if (changed) {
          if (!is_identity(opt_identity, get_value, old_kv)) {
            cl.send(old_kv, dest);
          }
          break;
        } else continue;
      }
    }
  }

  void send_with_tid(const arg_type& a, rank_type dest, int /*tid*/) {
    this->send(a, dest); // FIXME
  }

private:
  size_t get_bucket(const key_type& k) const {
    size_t h = hash_key(k);
    h ^= (h / this->get_size());
    h ^= (h / this->get_size() / this->get_size());
    h %= this->get_size();
    return h;
  }

  key_type get_dummy_value(size_t h) {
    for(size_t i = 0; ; ++i) {
      key_type x = dummy_value<key_type>()(i);
      if(get_bucket(x) != h)
	return x;
    }
    return dummy_value<key_type>()(0); // We should never be here
  }
};

template <typename CoalescingLayer, typename Combine, typename BufferSorter, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
class locking_cache_binop_reduction;

template <typename CoalescingLayerGen, typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
struct locking_cache_binop_reduction_gen {
  CoalescingLayerGen clg;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  unsigned int lg_size;

  public:
  locking_cache_binop_reduction_gen(const CoalescingLayerGen& clg, const Combine& combine, const OptIdentity& opt_identity, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval, unsigned int lg_size): clg(clg), combine(combine), opt_identity(opt_identity), get_key(get_key), get_value(get_value), make_keyval(make_keyval), lg_size(lg_size) {}

  template <typename Arg, typename Handler, typename BufferSorter = amplusplus::DummyBufferSorter<Arg> >
  struct inner {
    typedef locking_cache_binop_reduction<typename CoalescingLayerGen::template inner<Arg, Handler>::type, Combine, BufferSorter, OptIdentity, GetKey, GetValue, MakeKeyval> type;
  };
};

template <typename CoalescingLayer, typename Combine, typename BufferSorter, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval>
class locking_cache_binop_reduction {
  public:
  typedef typename message_type_traits<CoalescingLayer>::arg_type arg_type;
  typedef typename message_type_traits<CoalescingLayer>::handler_type handler_type;

  private:
  typedef typename boost::result_of<GetKey(arg_type)>::type key_type;
  typedef typename boost::result_of<GetValue(arg_type)>::type value_type;

  typedef transport::rank_type rank_type;

  CoalescingLayer cl;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  boost::hash<key_type> hash_key;
  unsigned int lg_size;
  // boost::scoped_array<amplusplus::detail::atomic<int /* bool */> > spinlocks; // atomic<bool> broken on BG/P
  // boost::scoped_array<arg_type> values;
  boost::scoped_array<std::pair<amplusplus::detail::atomic<int>, arg_type> > locks_and_values;
  amplusplus::detail::atomic<int> any_filled_entries;

  public:
  detail::hit_rate_counter counters;

  template <typename CoalescingLayerGen>
  explicit locking_cache_binop_reduction(const locking_cache_binop_reduction_gen<CoalescingLayerGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval>& gen, amplusplus::transport trans, valid_rank_set possible_dests, valid_rank_set possible_sources, const BufferSorter& bufstr)
    : cl(gen.clg, trans, possible_dests, possible_sources, bufstr),
      combine(gen.combine),
      opt_identity(gen.opt_identity),
      get_key(gen.get_key),
      get_value(gen.get_value),
      make_keyval(gen.make_keyval),
      hash_key(),
      lg_size(gen.lg_size),
      // spinlocks(new amplusplus::detail::atomic<int>[(size_t(1) << this->lg_size) * cl.get_transport().size()]),
      // values(new arg_type[(size_t(1) << this->lg_size) * cl.get_transport().size()]),
      locks_and_values(new std::pair<amplusplus::detail::atomic<int>, arg_type>[(size_t(1) << this->lg_size) * trans.size()]),
      any_filled_entries(0),
      counters(trans.get_nthreads())
  {
    this->clear();
  }

  size_t get_size() const {return size_t(1) << lg_size;}

  void clear() { // Clear without sending data
    counters.clear();
    size_t nranks = cl.get_transport().size();
    for (size_t i = 0; i < this->get_size() * nranks; ++i) {
      // spinlocks[i].store(false);
      locks_and_values[i].first.store(false);
      // values[i] = make_keyval(get_dummy_value(i % this->get_size()), value_type());
      locks_and_values[i].second = make_keyval(get_dummy_value(i % this->get_size()), value_type());
    }
    any_filled_entries.store(0);
  }

  static void lock_spinlock(amplusplus::detail::atomic<int>& l) {
    while (true) {
      int v = 0;
      if (l.compare_exchange_strong(v, 1)) return;
      BOOST_ASSERT (v == 1);
      amplusplus::detail::do_pause();
    }
  }

  static void unlock_spinlock(amplusplus::detail::atomic<int>& l) {
#ifdef NDEBUG
    l.store(0);
#else
    int old = l.exchange(0);
    BOOST_ASSERT (old == 1);
#endif
  }

  bool flush(rank_type dest) { // Clear while sending any valid data
    bool did_anything = false;
    for (size_t h = 0; h < this->get_size(); ++h) {
      size_t index = dest * this->get_size() + h;
      // if (spinlocks[index].load()) abort(); // Only safe for single-threaded use
      // lock_spinlock(spinlocks[index]);
      lock_spinlock(locks_and_values[index].first);
      // arg_type old_kv = values[index];
      arg_type old_kv = locks_and_values[index].second;
      // values[index] = make_keyval(get_dummy_value(h), value_type());
      locks_and_values[index].second = make_keyval(get_dummy_value(h), value_type());
      // unlock_spinlock(spinlocks[index]);
      unlock_spinlock(locks_and_values[index].first);

      key_type old_key = get_key(old_kv);
      size_t old_h = get_bucket(old_key);
      if (old_h == h && !is_identity(opt_identity, get_value, old_kv)) {
        cl.send(old_kv, dest);
        did_anything = true;
      }
    }
    any_filled_entries.store(0);
    return did_anything;
  }

  void send(const arg_type& a, rank_type dest) {
    key_type a_key = get_key(a);
    value_type a_value = get_value(a);
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.test(amplusplus::detail::get_thread_id());
#endif
    size_t h = get_bucket(hash_key(a_key));
    size_t index = dest * this->get_size() + h;
    // if (spinlocks[index].load()) abort(); // Only safe for single-threaded use
    // lock_spinlock(spinlocks[index]);
    lock_spinlock(locks_and_values[index].first);
    // arg_type old_kv = values[index];
    arg_type old_kv = locks_and_values[index].second;
    key_type old_key = get_key(old_kv);
    value_type old_value = get_value(old_kv);
    size_t old_h = get_bucket(old_key);
    if (a_key == old_key) {
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
      counters.hit(amplusplus::detail::get_thread_id());
#endif
      value_type new_value = combine(old_value, a_value);
      // values[index] = make_keyval(a_key, new_value);
      locks_and_values[index].second = make_keyval(a_key, new_value);
      // unlock_spinlock(spinlocks[index]);
      unlock_spinlock(locks_and_values[index].first);
    } else {
      // values[index] = a;
      locks_and_values[index].second = a;
      // unlock_spinlock(spinlocks[index]);
      unlock_spinlock(locks_and_values[index].first);
      if (old_h == h && !is_identity(opt_identity, get_value, old_kv)) {
        cl.send(old_kv, dest);
      } else {
        if (any_filled_entries.exchange(1) == 0) {
          cl.message_being_built(dest);
	  // This flush task may need to be fixed. We need flush tasks to be added once and then check if they have any work to do.
          cl.get_transport().add_flush_object(boost::bind(&locking_cache_binop_reduction::flush, this, dest));
        }
      }
    }
  }

  void send_with_tid(const arg_type& a, rank_type dest, int /*tid*/) {
    this->send(a, dest); // FIXME
  }

  const handler_type& get_handler() {return cl.get_handler();}
  void set_handler(const handler_type& h) {cl.set_handler(h);}
  transport get_transport() const {return cl.get_transport();}

  private:
  size_t get_bucket(const key_type& k) const {
    size_t h = hash_key(k);
    h ^= (h / this->get_size());
    h ^= (h / this->get_size() / this->get_size());
    h %= this->get_size();
    return h;
  }

  key_type get_dummy_value(size_t h) {
    for(size_t i = 0; ; ++i) {
      key_type x = dummy_value<key_type>()(i);
      if(get_bucket(x) != h)
	return x;
    }
    return dummy_value<key_type>()(0); // We should never be here
  }
};

template <typename CoalescingLayer, typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
class per_thread_cache_binop_reduction;

template <typename CoalescingLayerGen, typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
class per_thread_cache_binop_reduction_gen {
  CoalescingLayerGen clg;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  unsigned int lg_size;

  public:
  per_thread_cache_binop_reduction_gen(const CoalescingLayerGen& clg, const Combine& combine, const OptIdentity& opt_identity, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval, unsigned int lg_size): clg(clg), combine(combine), opt_identity(opt_identity), get_key(get_key), get_value(get_value), make_keyval(make_keyval), lg_size(lg_size) {}

  template <typename Arg, typename Handler>
  struct inner {
    typedef per_thread_cache_binop_reduction<CoalescingLayerGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval> type;
  };
};

template <typename CoalescingLayer, typename Combine, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval>
class per_thread_cache_binop_reduction {
  // BOOST_STATIC_WARNING((boost::is_same<CoalescingLayer, void>::value && !"per_thread_cache_binop_reduction is risky to use; use idempotent_combination instead of combination for operators such as min to get a write-through cache"));

  public:
  typedef typename message_type_traits<CoalescingLayer>::arg_type arg_type;
  typedef typename message_type_traits<CoalescingLayer>::handler_type handler_type;

  private:
  typedef typename boost::result_of<GetKey(arg_type)>::type key_type;
  typedef typename boost::result_of<GetValue(arg_type)>::type value_type;

  typedef transport::rank_type rank_type;

  CoalescingLayer cl;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  boost::hash<key_type> hash_key;
  unsigned int lg_size;
  valid_rank_set possible_dests;
  rank_type num_possible_dests;
  unsigned int num_threads;
  rank_type num_ranks;
  boost::scoped_array<arg_type> values;
  boost::scoped_array<arg_type*> start_ptrs;

  public:
  detail::hit_rate_counter counters;

  template <typename CoalescingLayerGen>
  explicit per_thread_cache_binop_reduction(const per_thread_cache_binop_reduction_gen<CoalescingLayerGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval>& gen, amplusplus::transport trans, valid_rank_set possible_dests = valid_rank_set(), valid_rank_set possible_sources = valid_rank_set())
    : cl(gen.clg, trans, possible_dests, possible_sources),
      combine(gen.combine),
      opt_identity(gen.opt_identity),
      get_key(gen.get_key),
      get_value(gen.get_value),
      make_keyval(gen.make_keyval),
      hash_key(),
      lg_size(gen.lg_size),
      possible_dests(possible_dests),
      num_possible_dests(possible_dests->count()),
      num_threads(trans.get_nthreads()),
      num_ranks(trans.size()),
      values(new arg_type[this->num_threads * this->num_possible_dests * (size_t(1) << this->lg_size)]),
      start_ptrs(new arg_type*[this->num_threads * this->num_ranks]),
      counters(trans.get_nthreads())
  {
    std::fill(&start_ptrs[0], &start_ptrs[num_threads * num_ranks],
              (arg_type*)0);
    for (size_t tid = 0; tid < num_threads; ++tid) {
      for (rank_type i = 0; i < num_possible_dests; ++i) {
        BOOST_ASSERT (possible_dests->rank_from_index(i) < num_ranks);
        start_ptrs[tid * num_ranks + possible_dests->rank_from_index(i)]
          = &values[(tid * num_possible_dests + i) * (size_t(1) << this->lg_size)];
      }
    }
    this->clear();
  }

  size_t get_size() const {return size_t(1) << lg_size;}

  void clear() { // Clear without sending data
    counters.clear();
    for (size_t i = 0; i < num_threads * num_possible_dests * this->get_size(); ++i) {
      values[i] = make_keyval(get_dummy_value(i % this->get_size()), value_type());
    }
  }

  bool flush(rank_type dest) { // Clear while sending any valid data
    bool did_anything = false;
    const int thread_id = amplusplus::detail::get_thread_id();
    BOOST_ASSERT (thread_id >= 0 && thread_id < num_threads);
    BOOST_ASSERT (dest < num_ranks);
    BOOST_ASSERT (start_ptrs[(size_t)thread_id * num_ranks + dest]);
    BOOST_ASSERT (start_ptrs[(size_t)thread_id * num_ranks + dest] >= values.get() + (size_t)thread_id * num_possible_dests * this->get_size());
    BOOST_ASSERT (start_ptrs[(size_t)thread_id * num_ranks + dest] < values.get() + (size_t)(thread_id + 1) * num_possible_dests * this->get_size());
    for (size_t h = 0; h < this->get_size(); ++h) {
      arg_type* ptr = start_ptrs[(size_t)thread_id * num_ranks + dest] + h;
      arg_type old_kv = *ptr;
      *ptr = make_keyval(get_dummy_value(h), value_type());

      key_type old_key = get_key(old_kv);
      size_t old_h = get_bucket(old_key);
      if (old_h == h && !is_identity(opt_identity, get_value, old_kv)) {
        cl.send_with_tid(old_kv, dest, thread_id);
        did_anything = true;
      }
    }
    return did_anything;
  }

  void send_with_tid(const arg_type& a, rank_type dest, int thread_id) {
    BOOST_STATIC_ASSERT (boost::is_void<arg_type>::value && !"Fix me up for new flush framework");
    key_type a_key = get_key(a);
    value_type a_value = get_value(a);
    counters.test(thread_id);
    BOOST_ASSERT (thread_id >= 0 && thread_id < num_threads);
    BOOST_ASSERT (dest < num_ranks);
    BOOST_ASSERT (start_ptrs[thread_id * num_ranks + dest]);
    BOOST_ASSERT (start_ptrs[(size_t)thread_id * num_ranks + dest] >= values.get() + (size_t)thread_id * num_possible_dests * this->get_size());
    BOOST_ASSERT (start_ptrs[(size_t)thread_id * num_ranks + dest] < values.get() + (size_t)(thread_id + 1) * num_possible_dests * this->get_size());
    size_t h = get_bucket(a_key);
    BOOST_ASSERT (h < this->get_size());
    arg_type* ptr = start_ptrs[(size_t)thread_id * num_ranks + dest] + h;
    arg_type old_kv = *ptr;
    key_type old_key = get_key(old_kv);
    value_type old_value = get_value(old_kv);
    if (a_key == old_key) {
      counters.hit(thread_id);
      value_type new_value = combine(old_value, a_value);
      *ptr = make_keyval(a_key, new_value);
    } else {
      *ptr = a;
      size_t old_h = get_bucket(old_key); // Check for dummy element
      if (old_h == h && !is_identity(opt_identity, get_value, old_kv)) {
        cl.send_with_tid(old_kv, dest, thread_id);
      }
    }
  }

  void send(const arg_type& a, rank_type dest) {
    this->send_with_tid(a, dest, amplusplus::detail::get_thread_id());
  }

  const handler_type& get_handler() {return cl.get_handler();}
  void set_handler(const handler_type& h) {cl.set_handler(h);}
  transport get_transport() const {return cl.get_transport();}

  private:
  size_t get_bucket(const key_type& k) const {
    size_t h = hash_key(k);
    h ^= (h / this->get_size());
    h ^= (h / this->get_size() / this->get_size());
    h %= this->get_size();
    return h;
  }

    key_type get_dummy_value(size_t h) {
    for(size_t i = 0; ; ++i) {
      key_type x = dummy_value<key_type>()(i);
      if(get_bucket(x) != h)
	return x;
    }
    return dummy_value<key_type>()(0); // We should never be here
  }
};

template <typename CoalescingLayerGen, typename Arg, typename Handler, typename Combine, typename BufferSorter = DummyBufferSorter<Arg>, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
class per_thread_wt_cache_binop_reduction;

template <typename CoalescingLayerGen, typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
struct per_thread_wt_cache_binop_reduction_gen {
  CoalescingLayerGen clg;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  unsigned int lg_size;

  public:
  per_thread_wt_cache_binop_reduction_gen(const CoalescingLayerGen& clg, const Combine& combine, const OptIdentity& opt_identity, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval, unsigned int lg_size): clg(clg), combine(combine), opt_identity(opt_identity), get_key(get_key), get_value(get_value), make_keyval(make_keyval), lg_size(lg_size) {}

  template <typename Arg, typename Handler, typename BufferSorter = amplusplus::DummyBufferSorter<Arg> >
  struct inner {
    typedef per_thread_wt_cache_binop_reduction<typename CoalescingLayerGen::template inner<Arg, Handler, BufferSorter>::type, Arg, Handler, Combine, BufferSorter, OptIdentity, GetKey, GetValue, MakeKeyval> type;
  };
};

template <typename CoalescingLayer, typename Arg, typename Handler, typename Combine, typename BufferSorter, typename OptIdentity, typename GetKey, typename GetValue, typename MakeKeyval>
class per_thread_wt_cache_binop_reduction {
  public:
  typedef Arg arg_type;
  typedef Handler handler_type;

  private:
  typedef typename boost::result_of<GetKey(arg_type)>::type key_type;
  typedef typename boost::result_of<GetValue(arg_type)>::type value_type;

  typedef transport::rank_type rank_type;

  CoalescingLayer cl;
  Combine combine;
  OptIdentity opt_identity;
  GetKey get_key;
  GetValue get_value;
  MakeKeyval make_keyval;
  boost::hash<key_type> hash_key;
  unsigned int lg_size;
  valid_rank_set possible_dests;
  rank_type num_possible_dests;
  unsigned int num_threads;
  rank_type num_ranks;
  boost::scoped_array<arg_type> values;
  boost::scoped_array<arg_type*> start_ptrs;

  public:
  detail::hit_rate_counter counters;

  template <typename CoalescingLayerGen>
  explicit per_thread_wt_cache_binop_reduction(const per_thread_wt_cache_binop_reduction_gen<CoalescingLayerGen, Combine, OptIdentity, GetKey, GetValue, MakeKeyval>& gen, amplusplus::transport trans, valid_rank_set possible_dests, valid_rank_set possible_sources, BufferSorter bufsrter)
    : cl(gen.clg, trans, possible_dests, possible_sources, bufsrter),
      combine(gen.combine),
      opt_identity(gen.opt_identity),
      get_key(gen.get_key),
      get_value(gen.get_value),
      make_keyval(gen.make_keyval),
      hash_key(),
      lg_size(gen.lg_size),
      possible_dests(possible_dests),
      num_possible_dests(possible_dests->count()),
      num_threads(cl.get_transport().get_nthreads()),
      num_ranks(cl.get_transport().size()),
      values(new arg_type[this->num_threads * this->num_possible_dests * (size_t(1) << this->lg_size)]),
      start_ptrs(new arg_type*[this->num_threads * this->num_ranks]),
      counters(trans.get_nthreads())
  {
    std::fill(&start_ptrs[0], &start_ptrs[num_threads * num_ranks],
              (arg_type*)0);
    for (size_t tid = 0; tid < num_threads; ++tid) {
      for (rank_type i = 0; i < num_possible_dests; ++i) {
        BOOST_ASSERT (possible_dests->rank_from_index(i) < num_ranks);
        start_ptrs[tid * num_ranks + possible_dests->rank_from_index(i)]
          = &values[(tid * num_possible_dests + i) * (size_t(1) << this->lg_size)];
      }
    }
    this->clear();
  }

  size_t get_size() const {return size_t(1) << lg_size;}

  void clear() { // Clear without sending data
    counters.clear();
    for (size_t i = 0; i < num_threads * num_possible_dests * this->get_size(); ++i) {
      values[i] = make_keyval(get_dummy_value(i % this->get_size()), value_type());
      BOOST_ASSERT (get_bucket(get_key(values[i])) != i % this->get_size());
    }
  }

  bool flush(rank_type dest) { // Clear while sending any valid data
    return false;
  }

  void send_with_tid(const arg_type& a, rank_type dest, const int thread_id) {
    key_type a_key = get_key(a);
    value_type a_value = get_value(a);
    counters.test(thread_id);
    BOOST_ASSERT (thread_id >= 0 && thread_id < (int)num_threads);
    BOOST_ASSERT (dest < num_ranks);
    BOOST_ASSERT (start_ptrs[thread_id * num_ranks + dest]);
    size_t h = get_bucket(a_key);
    BOOST_ASSERT (h < this->get_size());
    arg_type* ptr = start_ptrs[thread_id * num_ranks + dest] + h;
    arg_type old_kv = *ptr;
    key_type old_key = get_key(old_kv);
    value_type old_value = get_value(old_kv);
    value_type new_value = (a_key == old_key) ? combine(old_value, a_value) : a_value;
    if (a_key == old_key) {
      if (new_value != old_value) {
        *ptr = make_keyval(a_key, new_value);
        cl.send_with_tid(*ptr, dest, thread_id);
      } else {
	counters.hit(thread_id);
      }
    } else {
      *ptr = a;
      cl.send_with_tid(a, dest, thread_id);
    }
  }

  void send(const arg_type& a, rank_type dest) {
    this->send_with_tid(a, dest, amplusplus::detail::get_thread_id());
  }

  const handler_type& get_handler() {return cl.get_handler();}
  void set_handler(const handler_type& h) {cl.set_handler(h);}
  transport get_transport() const {return cl.get_transport();}

  private:
  size_t get_bucket(const key_type& k) const {
    size_t h = hash_key(k);
    h ^= (h / this->get_size());
    h ^= (h / this->get_size() / this->get_size());
    h %= this->get_size();
    return h;
  }

  key_type get_dummy_value(size_t h) {
    for(size_t i = 0; ; ++i) {
      key_type x = dummy_value<key_type>()(i);
      if(get_bucket(x) != h)
	return x;
    }
    return dummy_value<key_type>()(0); // We should never be here
  }
};

template <typename CoalescingLayerGen, typename ReductionGen, typename Policy>
class reduction_message_type_gen {
  typedef transport::rank_type rank_type;

  CoalescingLayerGen coalescing_layer_gen;
  ReductionGen reduction_gen;
  Policy policy;

  public:
  explicit reduction_message_type_gen(const CoalescingLayerGen& coalescing_layer_gen, const ReductionGen& reduction_gen, const Policy& policy = Policy()): coalescing_layer_gen(coalescing_layer_gen), reduction_gen(reduction_gen), policy(policy) {}

  template <typename Arg, typename Handler>
  class inner {
    public:
    typedef typename CoalescingLayerGen::template inner<Arg, Handler>::type coalescing_type;
    typedef typename CoalescingLayerGen::template inner<Arg, Handler>::gen_type coalescing_gen_type;
    typedef typename ReductionGen::template inner<coalescing_type, coalescing_gen_type, Policy>::type reduction_type;
    typedef typename ReductionGen::template inner<coalescing_type, coalescing_gen_type, Policy>::gen_type gen_type;
    typedef reduction_type type;
  };

  template <typename Arg, typename Handler>
  typename inner<Arg, Handler>::gen_type
  gen(transport trans, const valid_rank_set& possible_dests, const valid_rank_set& possible_sources) const {
    typename inner<Arg, Handler>::coalescing_gen_type cl(coalescing_layer_gen.template gen<Arg, Handler>(trans, possible_dests, possible_sources));
    return reduction_gen.template gen<typename inner<Arg, Handler>::coalescing_type, typename inner<Arg, Handler>::coalescing_gen_type, Policy>(cl, policy);
  }
};

#if 0
// template <typename Policy>
// class policy_from_binop {
//   Policy policy;

//   public:
//   explicit policy_from_binop(const Policy& policy = Policy()): policy(policy) {}

//   typedef typename Policy::pair_type value_type;

//   value_type combine(const value_type& a, const value_type& b) const {
//     if (policy.get_key(a) != policy.get_key(b)) {
//       return b;
//     } else {
//       typename Policy::value_type result = policy.op(policy.get_value(a), policy.get_value(b));
//       combined = a;
//       policy.get_value(combined) = result;
//     }
//   }

//   size_t hash(const value_type& a) const {
//     return policy.hash_key(policy.get_key(a));
//   }

//   value_type dummy_value(size_t h) const {
//     return policy.make_pair(policy.dummy_key(h), typename Policy::value_type());
//   }
// };
#endif

}

#endif // AMPLUSPLUS_REDUCTIONS_HPP
