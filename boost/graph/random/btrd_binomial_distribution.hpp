// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine

#ifndef BOOST_BTRD_BINOMIAL_HPP
#define BOOST_BTRD_BINOMIAL_HPP

#include <math.h>
#include <boost/random/geometric_distribution.hpp>

static double f_c(unsigned int k) {
  switch (k) {
    case 0: return .08106146679532726;
    case 1: return .04134069595540929;
    case 2: return .02767792568499834;
    case 3: return .02079067210376509;
    case 4: return .01664469118982119;
    case 5: return .01387612882307075;
    case 6: return .01189670994589177;
    case 7: return .01041126526197209;
    case 8: return .009255462182712733;
    case 9: return .008330563433362871;
    default: {
      double recip_k_plus_1 = 1. / (k + 1);
      double one_over_12 = 1./12;
      double one_over_360 = 1./360;
      double one_over_1260 = 1./1260;
      return (one_over_12 - (one_over_360 - one_over_1260 * recip_k_plus_1 * recip_k_plus_1) * recip_k_plus_1 * recip_k_plus_1) * recip_k_plus_1;
    }
  }
}

namespace boost { namespace graph { namespace random {

  template <typename Gen01>
  size_t btrd_binomial_distribution(size_t n_orig, double p, Gen01& gen) {
    // BTRD algorithm from pages 6--7 of "The Generation of Binomial Random
    // Variates" (Wolfgang Hoermann) --
    // http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.8407
    if (p > .5) return n_orig - btrd_binomial_distribution(n_orig, 1. - p, gen);
    if (n_orig * p < 10) {
      // First waiting time algorithm from page 525 of
      // http://cg.scs.carleton.ca/~luc/chapter_ten.pdf
      int x = 0;
      unsigned int sum = 0;
      do {
        sum += boost::geometric_distribution<>(1. - p)(gen);
        ++x;
      } while (sum <= n_orig);
      return x - 1;
    }
    if (n_orig > 1000000000) {
      return btrd_binomial_distribution(1000000000, p, gen) +
             btrd_binomial_distribution(n_orig - 1000000000, p, gen);
    }
    int n = (int)n_orig;
    // Setup
    int m = (int)floor((n + 1) * p);
    double r = p / (1. - p);
    double nr = (n + 1) * r;
    double npq = n * p * (1. - p);
    double b = 1.15 + 2.53 * sqrt(npq);
    double a = -.0873 + .0248 * b + .01 * p;
    double c = n * p + .5;
    double alpha = (2.83 + 5.1 / b) * sqrt(npq);
    double v_r = .92 - 4.2 / b;
    double u_rv_r = .86 * v_r;
    while (true) {
      // 1
      double v = gen();
      double u;
      if (v <= u_rv_r) {
        u = v / v_r - .43;
        return (int)floor((2 * a / (.5 + fabs(u)) + b) * u + c);
      }
      // 2
      if (v >= v_r) {
        u = gen() - .5;
      } else {
        u = v / v_r - .93;
        u = copysign(1., u) * .5 - u;
        v = v_r * gen();
      }
      // 3.0
      double us = .5 - fabs(u);
      int k = (int)floor((2 * a / us + b) * u + c);
      if (k < 0 || k > n) continue;
      v *= alpha / (a / (us * us) + b);
      int km = abs(k - m);
      if (km > 15) {
        // 3.2
        v = log(v);
        double rho = (km / npq) * (((km / 3 + .625) * km + 1. / 6) / npq + .5);
        double t = -km * km / (2 * npq);
        if (v < t - rho) return k;
        if (v > t + rho) continue;
        // 3.3
        int nm = n - m + 1;
        double h = (m + .5) * log((m + 1) / (r * nm)) + f_c(m) + f_c(n - m);
        // 3.4
        int nk = n - k + 1;
        double threshold = h +
                           (n + 1) * log((double)nm / nk) +
                           (k + .5) * log(nk * r / (k + 1)) -
                           f_c(k) -
                           f_c(n - k);
        if (v <= threshold) return k;
      } else {
        // 3.1
        double f = 1.;
        if (m < k) {
          for (int i = m; i != k; ++i) f *= nr / i - r;
        } else if (m > k) {
          for (int i = k; i != m; ++i) v *= nr / i + r;
        }
        if (v <= f) return k;
      }
    }
  }

} } } // end namespace boost::graph::random 

#endif // BOOST_BTRD_BINOMIAL_HPP
