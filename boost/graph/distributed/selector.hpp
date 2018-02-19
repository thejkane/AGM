// Copyright (C) 2006 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_DISTRIBUTED_SELECTOR_HPP
#define BOOST_GRAPH_DISTRIBUTED_SELECTOR_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/version.hpp>

#if BOOST_VERSION >= 105300
#include <boost/graph/detail/is_distributed_selector.hpp>
#endif

namespace boost { 

  /* The default local selector for a distributedS selector. */
  struct defaultS {};

  /**
   * Selector that specifies that the graph should be distributed
   * among different processes.
   */
  template<typename LocalS = defaultS, typename DistributionS = defaultS>
  struct distributedS 
  {
    typedef LocalS local_selector;
    typedef DistributionS distribution;
  };

#if BOOST_VERSION >= 105300
  namespace detail {
    template <typename LocalS, typename DistributionS>
    struct is_distributed_selector<distributedS<LocalS, DistributionS> >: boost::mpl::true_ {};
  }
#endif
}

#endif // BOOST_GRAPH_DISTRIBUTED_SELECTOR_HPP
