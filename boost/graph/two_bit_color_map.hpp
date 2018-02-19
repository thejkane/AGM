// Copyright (C) 2005-2012 The Trustees of Indiana University.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Jeremiah Willcock
//           Douglas Gregor
//           Nicholas Edmonds
//           Andrew Lumsdaine

// Two bit per color property map

#ifndef BOOST_TWO_BIT_COLOR_MAP_HPP
#define BOOST_TWO_BIT_COLOR_MAP_HPP

#include <boost/property_map/property_map.hpp>
#include <boost/shared_array.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <algorithm>

namespace boost {

enum two_bit_color_type { 
  two_bit_white = 0, 
  two_bit_gray  = 1, 
  two_bit_green = 2, 
  two_bit_black = 3 
};

template <>
struct color_traits<two_bit_color_type>
{
  static two_bit_color_type white() { return two_bit_white; }
  static two_bit_color_type gray()  { return two_bit_gray; }
  static two_bit_color_type green() { return two_bit_green; }
  static two_bit_color_type black() { return two_bit_black; }
};


template<typename IndexMap = identity_property_map>
struct two_bit_color_map 
{
  typedef unsigned int word_type;

  std::size_t n;
  IndexMap index;
  static const std::size_t BITS_PER_WORD = sizeof(word_type) * CHAR_BIT;  // if CHAR_BIT is odd we just ignore the highest bit
  shared_array<word_type> data;

  typedef typename property_traits<IndexMap>::key_type key_type;
  typedef two_bit_color_type value_type;
  typedef void reference;
  typedef read_write_property_map_tag category;

  explicit two_bit_color_map(std::size_t n, const IndexMap& index = IndexMap())
    : n(n), index(index), data(new word_type[(n + (BITS_PER_WORD / 2 - 1)) / (BITS_PER_WORD / 2)])
  {
    // Fill to white
    std::fill(data.get(), data.get() + (n + (BITS_PER_WORD / 2 - 1)) / (BITS_PER_WORD / 2), 0);
  }
};

template<typename IndexMap>
inline two_bit_color_type
get(const two_bit_color_map<IndexMap>& pm, 
    typename two_bit_color_map<IndexMap>::key_type key) 
{
  typename property_traits<IndexMap>::value_type i = get(pm.index, key);
  assert ((std::size_t)i < pm.n);
  return two_bit_color_type((pm.data.get()[i / (pm.BITS_PER_WORD / 2)] >> ((i % (pm.BITS_PER_WORD / 2)) * 2)) & 3);
}

template<typename IndexMap>
inline void
put(const two_bit_color_map<IndexMap>& pm, 
    typename two_bit_color_map<IndexMap>::key_type key,
    two_bit_color_type value)
{
  typedef typename two_bit_color_map<IndexMap>::word_type word_type; 

  typename property_traits<IndexMap>::value_type i = get(pm.index, key);
  assert ((std::size_t)i < pm.n);
  //assert (value >= 0 && value < 4);
  std::size_t word_num = i / (pm.BITS_PER_WORD / 2);
  std::size_t bit_position = ((i % (pm.BITS_PER_WORD / 2)) * 2);
  word_type* word_ptr = &pm.data.get()[word_num];
  *word_ptr = (*word_ptr & ~((word_type)3 << bit_position))
            | ((word_type)value << bit_position);
}

template<typename IndexMap>
inline bool
exchange(const two_bit_color_map<IndexMap>& pm, 
	 typename two_bit_color_map<IndexMap>::key_type key,
	 two_bit_color_type oldval, two_bit_color_type newval)
{
  using boost::parallel::val_compare_and_swap;

  typedef typename two_bit_color_map<IndexMap>::word_type word_type; 

  typename property_traits<IndexMap>::value_type i = get(pm.index, key);
  assert ((std::size_t)i < pm.n);
  // assert (newval >= 0 && newval < 4);
  std::size_t word_num = i / (pm.BITS_PER_WORD / 2);
  std::size_t bit_position = ((i % (pm.BITS_PER_WORD / 2)) * 2);
  word_type* word_ptr = &pm.data.get()[word_num];
  word_type old_word = *word_ptr;
  word_type new_word, last_old_word;

  while ((two_bit_color_type)((old_word >> bit_position) & 3) == oldval) {
    last_old_word = old_word;
    new_word = (old_word & ~((word_type)3 << bit_position)) | ((word_type)newval << bit_position);
    old_word = val_compare_and_swap(word_ptr, old_word, new_word);
    if (old_word == last_old_word)
      return true;
  }
  
  return false;
}

template<typename IndexMap>
inline two_bit_color_map<IndexMap>
make_two_bit_color_map(std::size_t n, const IndexMap& index_map)
{
  return two_bit_color_map<IndexMap>(n, index_map);
}

} // end namespace boost

#endif // BOOST_TWO_BIT_COLOR_MAP_HPP

#ifdef BOOST_GRAPH_USE_MPI
#  include <boost/graph/distributed/two_bit_color_map.hpp>
#endif
