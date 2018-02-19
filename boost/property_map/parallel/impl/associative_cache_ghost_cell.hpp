// Copyright (C) 2004-2009 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Nick Edmonds
//           Andrew Lumsdaine

// The placement of this #include probably looks very odd relative to
// the #ifndef/#define pair below. However, this placement is
// extremely important to allow the various property map headers to be
// included in any order.
#include <boost/property_map/property_map.hpp>

#ifndef BOOST_PARALLEL_IMPL_ASSOCIATIVE_CACHE_GHOST_CELL_HPP
#define BOOST_PARALLEL_IMPL_ASSOCIATIVE_CACHE_GHOST_CELL_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/graph/parallel/detail/untracked_pair.hpp>
#include <vector>
#include <algorithm>
#include <boost/iterator/filter_iterator.hpp>

namespace boost {
  namespace parallel {
    namespace detail {

template <typename Key, typename Value>
struct cache_entry {
  Key key;
  Value value;
  bool valid;
  cache_entry(): key(), value(), valid(false) {}
};

struct get_cache_entry_valid {
  template <typename Key, typename Value>
  bool operator()(const cache_entry<Key, Value>& e) const {
    return e.valid;
  }
};

template <typename Key, typename Value, int Associativity>
class associative_cache_ghost_cell_storage {
  typedef cache_entry<Key, Value> entry;
  std::vector<entry> data;

  public:
  /// Iterator into the ghost cells
  typedef boost::filter_iterator<get_cache_entry_valid, typename std::vector<entry>::iterator> iterator;

  static const Key& get_key(iterator i) {return i->key;}
  static Value& get_value(iterator i) {return i->value;}

  void get_cell(const Key& key, Value*& value_ptr, bool& found) {
    std::size_t key_n = boost::hash<Key>()(key);
    entry* bucket = &data[key_n % (data.size() / Associativity) * Associativity];
    for (size_t i = 0; i < Associativity; ++i) {
      if (bucket[i].valid && bucket[i].key == key) {
        found = true;
        if (i != 0) {
          // Move this entry to the front
          std::rotate(&bucket[0], &bucket[i], &bucket[Associativity]);
        }
        value_ptr = &bucket[0].value;
        return;
      }
    }
    found = false;
    value_ptr = NULL;
    return;
  }

  template <typename PMap>
  Value& insert_cell(const Key& key, const Value& value, PMap* pm) {
    std::size_t key_n = boost::hash<Key>()(key);
    entry* bucket = &data[key_n % (data.size() / Associativity) * Associativity];
    for (size_t i = 0; i < Associativity; ++i) {
      if (!bucket[i].valid) {
        bucket[i].key = key;
        bucket[i].value = value;
        bucket[i].valid = true;
        if (i != 0) {
          // Move this entry to the front
          std::rotate(&bucket[0], &bucket[i], &bucket[Associativity]);
        }
        return bucket[0].value;
      }
    }
    // Purge last (LRU) value
    if (Associativity != 1) {
      std::rotate(&bucket[0], &bucket[Associativity - 1], &bucket[Associativity]);
    }
    pm->prune_one_ghost_cell(bucket[0].key, bucket[0].value);
    bucket[0].key = key;
    bucket[0].value = value;
    bucket[0].valid = true;
    return bucket[0].value;
  }

  template <typename PMap>
  void prune(PMap* pm) {}

  template <typename PMap>
  void set_limit(size_t m, PMap* pm) {
    // Resizing flushes all entries
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i].valid) {
        pm->prune_one_ghost_cell(data[i].key, data[i].value);
        data[i].valid = false;
      }
    }
    // Size is rounded up to multiple of Associativity
    data.resize((m + Associativity - 1) / Associativity * Associativity);
  }

  void unset_limit() {
    assert (!"Cannot unset size limit");
  }

  void clear() {
    data.clear();
  }

  iterator begin() {
    return iterator(get_cache_entry_valid(), data.begin());
  }

  iterator end() {
    return iterator(get_cache_entry_valid(), data.end());
  }
};

    }
  }
}

#endif // BOOST_PARALLEL_IMPL_ASSOCIATIVE_CACHE_GHOST_CELL_HPP
