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

#ifndef BOOST_PARALLEL_IMPL_MULTI_INDEX_GHOST_CELL_HPP
#define BOOST_PARALLEL_IMPL_MULTI_INDEX_GHOST_CELL_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/graph/parallel/detail/untracked_pair.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/sequenced_index.hpp>

namespace boost {
  namespace parallel {
    namespace detail {

template <typename Key, typename Value>
class multi_index_ghost_cell_storage {
  /// The type of the ghost cells
  typedef multi_index::multi_index_container<
            std::pair<Key, Value>,
            multi_index::indexed_by<
              multi_index::sequenced<>,
              multi_index::hashed_unique<
                multi_index::member<std::pair<Key, Value>,
                                    Key,
                                    &std::pair<Key, Value>::first>
              >
            >
          > ghost_cells_type;

  public:
  /// Iterator into the ghost cells
  typedef typename ghost_cells_type::iterator iterator;

  static const Key& get_key(iterator i) {return i->first;}
  static Value& get_value(iterator i) {return const_cast<Value&>(i->second);}

  private:
  /// Key-based index into the ghost cells
  typedef typename ghost_cells_type::template nth_index<1>::type
    ghost_cells_key_index_type;

  /// Iterator into the ghost cells (by key)
  typedef typename ghost_cells_key_index_type::iterator key_iterator;

  ghost_cells_type data;
  size_t max_ghost_cells;

  public:
  void get_cell(const Key& key, Value*& value_ptr, bool& found) {
    // Index by key
    ghost_cells_key_index_type const& key_index 
      = data.template get<1>();

    // Search for the ghost cell by key, and project back to the sequence
    iterator ghost_cell 
      = data.template project<0>(key_index.find(key));
    if (ghost_cell == data.end()) {
      found = false;
      value_ptr = 0;
    } else {
      if (max_ghost_cells > 0) {
        // Put this cell at the beginning of the MRU list
        data.relocate(data.begin(), ghost_cell);
      }
      found = true;
      value_ptr = const_cast<Value*>(&ghost_cell->second);
    }
  }

  template <typename PMap>
  Value& insert_cell(const Key& key, const Value& value, PMap* pm) {
    // Create a ghost cell containing the new value
    iterator ghost_cell = data.push_front(std::make_pair(key, value)).first;

    // If we need to, prune the ghost cells
    if (max_ghost_cells > 0)
      this->prune(pm);

    return const_cast<Value&>(ghost_cell->second);
  }

  template <typename PMap>
  void prune(PMap* pm) {
    if (max_ghost_cells == 0)
      return;

    while (data.size() > max_ghost_cells) {
      // Evict the last ghost cell
      pm->prune_one_ghost_cell(data.back().first, data.back().second);

      // Actually remove the ghost cell
      data.pop_back();
    }
  }

  template <typename PMap>
  void set_limit(size_t m, PMap* pm) {
    if (m == 1)
      // It is not safe to have only 1 ghost cell; the cell() method
      // will fail.
      m = 2;

    max_ghost_cells = m;
    this->prune(pm);
  }

  void unset_limit() {
    max_ghost_cells = 0;
  }

  void clear() {
    data.clear();
  }

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

#endif // BOOST_PARALLEL_IMPL_MULTI_INDEX_GHOST_CELL_HPP
