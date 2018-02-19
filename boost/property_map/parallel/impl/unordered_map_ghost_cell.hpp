// Copyright (C) 2004-2009 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Andrew Lumsdaine

// The placement of this #include probably looks very odd relative to
// the #ifndef/#define pair below. However, this placement is
// extremely important to allow the various property map headers to be
// included in any order.
#include <boost/property_map/property_map.hpp>

#ifndef BOOST_PARALLEL_IMPL_UNORDERED_MAP_GHOST_CELL_HPP
#define BOOST_PARALLEL_IMPL_UNORDERED_MAP_GHOST_CELL_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/unordered_map.hpp>

namespace boost {
  namespace parallel {
    namespace detail {

template <typename Key, typename Value>
class unordered_map_ghost_cell_storage {
  // The type of the ghost cells
  typedef unordered_map<Key, Value> ghost_cells_type;

public:

  // Iterator into the ghost cells
  typedef typename ghost_cells_type::iterator iterator;

  static const Key& get_key(iterator i) {return i->first;}
  static Value& get_value(iterator i) {return const_cast<Value&>(i->second);}

private:

  ghost_cells_type data;

public: 

  void get_cell(const Key& key, Value*& value_ptr, bool& found) {
    iterator loc;

    if ((loc = data.find(key)) == data.end()) {
      found = false;
      value_ptr = 0;
    } else {
      found = true;
      value_ptr = &loc->second;
    }
  }

  template <typename PMap>
  Value& insert_cell(const Key& key, const Value& value, PMap* pm) {
    iterator loc = (data.insert(std::make_pair(key, value))).first;
    return loc->second;
  }

  template <typename PMap>
  void prune(PMap* pm) { assert(false); }

  template <typename PMap>
  void set_limit(size_t m, PMap* pm) { }

  void unset_limit() {}

  void clear () { data.clear(); }

  iterator begin() {
    return data.begin();
  }

  iterator end() {
    return data.end();
  }
};

    }
  }
}

#endif // BOOST_PARALLEL_IMPL_UNORDERED_MAP_GHOST_CELL_HPP
