// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#ifndef BOOST_NWAY_BERNOULLI_DISTRIBUTION_HPP
#define BOOST_NWAY_BERNOULLI_DISTRIBUTION_HPP

#include <boost/random/uniform_01.hpp>
#include <limits>

namespace boost { namespace graph { namespace random {

  template <int NCases, typename RAIter>
  struct find_bucket_impl {
    typedef typename std::iterator_traits<RAIter>::value_type value_type;
    static int go(value_type x, RAIter start) {
      RAIter middle = start + NCases / 2;
      if (x < *middle) {
        return find_bucket_impl<NCases / 2, RAIter>::go(x, start);
      } else {
        return NCases / 2 + find_bucket_impl<NCases - NCases / 2, RAIter>::go(x, middle);
      }
    }
  };

  template <typename RAIter>
  struct find_bucket_impl<1, RAIter> {
    typedef typename std::iterator_traits<RAIter>::value_type value_type;
    static int go(value_type /* x */, RAIter /* start */) {
      return 0;
    }
  };

  template <typename RAIter>
  struct find_bucket_impl<0, RAIter> {};

  // Find the last bucket whose value is <= x
  template <int NCases, typename RAIter>
  int find_bucket(typename std::iterator_traits<RAIter>::value_type x,
                  RAIter start) {
    return find_bucket_impl<NCases, RAIter>::go(x, start);
  }

  template <int NCases, typename RealType = double, typename RNGResult = RealType>
  class nway_bernoulli_distribution {
    RNGResult ranges[NCases + 1];

    public:
    nway_bernoulli_distribution() {}

    // Range of generator must be [min_, max_)
    template <typename RandomAccessIterator>
    explicit nway_bernoulli_distribution(RandomAccessIterator it, RNGResult min_, RNGResult max_) {
      set_probabilities(it, min_, max_);
    }

    template <typename RandomAccessIterator>
    void set_probabilities(RandomAccessIterator it, RNGResult min_, RNGResult max_) {
      RealType total = 0;
      for (size_t i = 0; i < NCases; ++i) {
	ranges[i] = RNGResult(min_ + (max_ - min_) * total);
	total += *it++;
      }
      ranges[NCases] = max_;
    }

    int compute_from_random_number(RNGResult x) const {
      return find_bucket<NCases>(x, ranges);
    }

    template <typename Gen>
    int operator()(Gen& gen) const {
      return compute_from_random_number(gen());
    }
  };

  template <int NCases, typename RealType, typename RandomAccessIterator, typename Gen>
  nway_bernoulli_distribution<NCases, RealType, typename Gen::result_type>
  make_nway_bernoulli_distribution(RandomAccessIterator it, const Gen& gen) {
    return nway_bernoulli_distribution<NCases, RealType, typename Gen::result_type>
             (it,
              (gen.min)(),
              (gen.max)() + std::numeric_limits<typename Gen::result_type>::is_integer);
  }

} } } // end namespace boost::graph::random

#endif // BOOST_NWAY_BERNOULLI_DISTRIBUTION_HPP
