// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

// The placement of this #include probably looks very odd relative to
// the #ifndef/#define pair below. However, this placement is
// extremely important to allow the various property map headers to be
// included in any order.
#include <boost/property_map/property_map.hpp>

#ifndef BOOST_PARALLEL_LOCAL_PROPERTY_MAP_HPP
#define BOOST_PARALLEL_LOCAL_PROPERTY_MAP_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <cassert>

namespace boost {
  /** Property map that accesses an underlying, local property map
   * using a subset of the global keys.
   */
  template<typename GlobalMap, typename StorageMap>
  class local_property_map
  {
    typedef typename property_traits<GlobalMap>::value_type owner_local_pair;

  public:
    typedef typename property_traits<StorageMap>::value_type value_type;
    typedef typename property_traits<GlobalMap>::key_type key_type;
    typedef typename property_traits<StorageMap>::reference  reference;
    typedef typename property_traits<StorageMap>::category   category;

    local_property_map() { }

    local_property_map(amplusplus::transport& trans, 
                       const GlobalMap& global, const StorageMap& storage)
      : transport_(trans), global_(global), storage(storage) { }

    // Cast away the const'ness of the transport since the const graph it comes
    // from is in sequential BGL and I don't want to pull in that whole library 
    // to modify it -- this is not safe, but it is expedient (NGE)
    local_property_map(const amplusplus::transport& trans, 
                       const GlobalMap& global, const StorageMap& storage)
      : transport_(const_cast<amplusplus::transport&>(trans)), global_(global), storage(storage) { }

    reference operator[](const key_type& key) 
    { 
      owner_local_pair p = get(global_, key);
      assert(p.first == transport_.rank());
      return storage[p.second]; 
    }

    GlobalMap& global() const { return global_; }
    StorageMap& base() const { return storage; }

    amplusplus::transport& transport()       { return transport_; }
    amplusplus::transport& transport() const { return transport_; }

  private:
    amplusplus::transport& transport_;
    mutable GlobalMap global_;
    mutable StorageMap storage;
  };

  template<typename GlobalMap, typename StorageMap>
  inline
  typename local_property_map<GlobalMap, StorageMap>::reference
  get(const local_property_map<GlobalMap, StorageMap>& pm, 
      typename local_property_map<GlobalMap, StorageMap>::key_type
        const & key)

  {
    typename property_traits<GlobalMap>::value_type p = get(pm.global(), key);
    return get(pm.base(), p.second);
  }

  template<typename GlobalMap, typename StorageMap>
  inline void
  put(const local_property_map<GlobalMap, StorageMap>& pm, 
      typename local_property_map<GlobalMap, StorageMap>
                 ::key_type const & key,
      typename local_property_map<GlobalMap, StorageMap>
                 ::value_type const& v)
  {
    typename property_traits<GlobalMap>::value_type p = get(pm.global(), key);
    assert(p.first == pm.transport().rank());
    put(pm.base(), p.second, v);
  }

  template<typename GlobalMap, typename StorageMap>
  inline bool
  exchange(const local_property_map<GlobalMap, StorageMap>& pm, 
	   typename local_property_map<GlobalMap, StorageMap>
	              ::key_type const & key,
	   typename local_property_map<GlobalMap, StorageMap>
	              ::value_type const& oldval,
	   typename local_property_map<GlobalMap, StorageMap>
                      ::value_type const& newval)
  {
    typename property_traits<GlobalMap>::value_type p = get(pm.global(), key);
    assert(p.first == pm.transport().rank());
    return exchange(pm.base(), p.second, oldval, newval);
  }

  template<typename GlobalMap, typename StorageMap>
  inline bool
  maybe_exchange(const local_property_map<GlobalMap, StorageMap>& pm, 
		 typename local_property_map<GlobalMap, StorageMap>
	                    ::key_type const & key,
		 typename local_property_map<GlobalMap, StorageMap>
	                    ::value_type const& oldval,
		 typename local_property_map<GlobalMap, StorageMap>
                            ::value_type const& newval)
  {
    typename property_traits<GlobalMap>::value_type p = get(pm.global(), key);
    assert(p.first == pm.transport().rank());
    return maybe_exchange(pm.base(), p.second, oldval, newval);
  }
} // end namespace boost
#endif // BOOST_PARALLEL_LOCAL_PROPERTY_MAP_HPP
