// Copyright (C) 2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine
#ifndef BOOST_PARALLEL_GRAPH_UTILITY_HPP
#define BOOST_PARALLEL_GRAPH_UTILITY_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <functional>

namespace boost { namespace parallel {

  template <typename T>
  struct identity {
    typedef T result_type;
    T operator() (const T& t) const { return t; }
  };

  //
  // A property map indexed by pairs which maps pair.first as the index to PM
  //
  struct map_of_1st {
    typedef amplusplus::transport::rank_type rank_type;
    
    template <typename Snd>
    friend rank_type get(const map_of_1st& m, const std::pair<rank_type, Snd>& p)
    { return p.first; }
  };
  
  map_of_1st make_map_of_1st() { return map_of_1st(); }

  //
  // A property map indexed by pairs which maps pair.first as the index to PM
  //
  template <typename PM>
  struct map_of_project1st {
    typedef typename boost::property_traits<PM>::value_type value_type;
    typedef typename boost::property_traits<PM>::key_type key_first_type;
    
    map_of_project1st(const PM& pm) : pm(pm) {}
    
    template <typename Snd>
    friend value_type get(const map_of_project1st& m,
			  const std::pair<key_first_type, Snd>& p)
    { return get(m.pm, p.first); }
    
  private:
    PM pm;
  };
  
  template <typename PM>
  map_of_project1st<PM> make_map_of_project1st(const PM& pm) {
    return map_of_project1st<PM>(pm);
  }

  //
  // A property map indexed by pairs which maps pair.first as the index to PM
  // then strips off the first element of the return type (needed for distributed PMap)
  // NGE: There has to be a more elegant way to do this probably with Boost.Lambda
  //
  template <typename PM>
  struct first_of_map_of_project1st {
    typedef typename boost::property_traits<PM>::value_type::first_type value_type;
    typedef typename boost::property_traits<PM>::key_type key_first_type;
    
    first_of_map_of_project1st(const PM& pm) : pm(pm) {}
    
    template <typename Snd>
    friend value_type get(const first_of_map_of_project1st& m,
			  const std::pair<key_first_type, Snd>& p)
    { return get(m.pm, p.first).first; }
    
  private:
    PM pm;
  };
  
  template <typename PM>
  first_of_map_of_project1st<PM> make_first_of_map_of_project1st(const PM& pm) {
    return first_of_map_of_project1st<PM>(pm);
  }

  //
  // A property map that returns the first element of the underlying property map
  // NGE: Ditto Boost.Lamda's likely usefulness here
  //
  template <typename PM>
  struct first_of_map {
    typedef typename boost::property_traits<PM>::value_type::first_type value_type;
    typedef typename boost::property_traits<PM>::key_type key_first_type;
    
    first_of_map(const PM& pm) : pm(pm) {}
    
    friend value_type get(const first_of_map& m,
			  const key_first_type& p)
    { return get(m.pm, p).first; }
    
  private:
    PM pm;
  };
  
  template <typename PM>
  first_of_map<PM> make_first_of_map(const PM& pm) {
    return first_of_map<PM>(pm);
  }

  //
  // A property map indexed by pairs which maps a given tuple element as the index to PM
  //  
  template <typename PM, int Idx>
  struct map_of_tuple_elem {
    typedef typename boost::property_traits<PM>::value_type value_type;
    typedef typename boost::property_traits<PM>::key_type key_first_type;
    
    map_of_tuple_elem(const PM& pm): pm(pm) {}
    
    template <typename Tuple>
    friend value_type get(const map_of_tuple_elem& m, const Tuple& p)
    { return get(m.pm, p.template get<Idx>()); }
    
  private:
    PM pm;
  };
  
  template <int Idx, typename PM>
  map_of_tuple_elem<PM, Idx> make_map_of_tuple_elem(const PM& pm) {
    return map_of_tuple_elem<PM, Idx>(pm);
  }

  //
  // A property map indexed by pairs which maps pair.first as the index to PM
  //
  template <int Idx>
  struct tuple_elem {
    typedef amplusplus::transport::rank_type rank_type;
    
    template <typename Tuple>
    friend rank_type get(const tuple_elem& m, const Tuple& p)
    { return p.template get<Idx>(); }
  };
  
  template <int Idx>
  tuple_elem<Idx> make_tuple_elem() { return tuple_elem<Idx>(); }

} } // boost::parallel

#endif // BOOST_PARALLEL_GRAPH_UTILITY_HPP
