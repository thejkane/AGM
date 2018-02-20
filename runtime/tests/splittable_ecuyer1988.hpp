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

#ifndef BOOST_SPLITTABLE_ECUYER1988_HPP
#define BOOST_SPLITTABLE_ECUYER1988_HPP

#include <boost/cstdint.hpp>
#include <boost/iterator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <utility>

namespace boost {
  namespace random {

    template <typename Base> // To get around incomplete class issues
    class generic_split_iterator: public boost::iterator_facade<generic_split_iterator<Base>, Base, std::input_iterator_tag, Base> {
      public:
      Base base;
      uint32_t new_multiplier1, new_multiplier2;

      generic_split_iterator(Base base, uint32_t new_multiplier1, uint32_t new_multiplier2)
	: base(base), new_multiplier1(new_multiplier1), new_multiplier2(new_multiplier2) {}

      generic_split_iterator(): base(), new_multiplier1(0), new_multiplier2(0) {}

      Base dereference() const {
	Base result = base;
	result.multiplier1 = new_multiplier1;
	result.multiplier2 = new_multiplier2;
	return result;
      }

      bool equal(const generic_split_iterator& it) const {
	return base.state1 == it.base.state1 && base.state2 == it.base.state2;
      }

      void increment() {
	base(); // Advance seeds
      }
    };

    class splittable_ecuyer1988 {
      uint32_t state1, state2, multiplier1, multiplier2;

      static const uint64_t modulus1 = 2147483563;
      static const uint64_t modulus2 = 2147483399;

      template <uint64_t Modulus>
      static uint32_t mul_mod(uint32_t a, uint32_t b) {
	return uint32_t((uint64_t(a) * b) % Modulus);
      }

      template <uint64_t Modulus>
      static uint32_t exp_mod(uint32_t a, uint32_t b) {
	if (b == 0) {
	  return 1;
	} else if (b == 1) {
	  return a;
	} else {
	  uint32_t temp = exp_mod<Modulus>(a, b>>1);
	  temp = mul_mod<Modulus>(temp, temp);
	  if (b % 2) temp = mul_mod<Modulus>(a, temp);
	  return temp;
	}
      }

      public:
      splittable_ecuyer1988(): state1(1), state2(1), multiplier1(40014), multiplier2(40692) {}
      splittable_ecuyer1988(uint32_t seed1, uint32_t seed2): state1(seed1), state2(seed2), multiplier1(40014), multiplier2(40692) {}

      void seed(uint32_t seed1, uint32_t seed2) {
	state1 = seed1;
	state2 = seed2;
      }

      void split(splittable_ecuyer1988& a, splittable_ecuyer1988& b) const {
	a = *this;
	a.multiplier1 = mul_mod<modulus1>(a.multiplier1, a.multiplier1);
	a.multiplier2 = mul_mod<modulus2>(a.multiplier2, a.multiplier2);
	b = *this;
	b();
	b.multiplier1 = a.multiplier1;
	b.multiplier2 = a.multiplier2;
      }

      typedef generic_split_iterator<splittable_ecuyer1988> split_iterator;
      typedef std::pair<split_iterator, split_iterator> split_iterator_pair;
      friend class generic_split_iterator<splittable_ecuyer1988>;

      split_iterator_pair split_n(uint32_t n) const { // This RNG is not usable after the split
	split_iterator begin(*this, exp_mod<modulus1>(multiplier1, n), exp_mod<modulus2>(multiplier2, n));
	split_iterator end(*this, begin.new_multiplier1, begin.new_multiplier2);
	end.base.state1 = mul_mod<modulus1>(end.base.state1, end.new_multiplier1);
	end.base.state2 = mul_mod<modulus2>(end.base.state2, end.new_multiplier2);
	return std::make_pair(begin, end);
      }

      split_iterator_pair split_off_n(uint32_t n) { // Allows this RNG to keep being used
	split_iterator_pair p = split_n(n + 1);
	*this = *p.first;
	++p.first;
	return p;
      }

      typedef uint32_t result_type;

      uint32_t operator()() {
	state1 = mul_mod<modulus1>(state1, multiplier1);
	state2 = mul_mod<modulus2>(state2, multiplier2);
	return (state1 + modulus1 - 1 - state2) % (modulus1 - 1);
      }

      BOOST_STATIC_CONSTANT(bool, has_fixed_range = true);
      BOOST_STATIC_CONSTANT(uint32_t, min_value = 0);
      BOOST_STATIC_CONSTANT(uint32_t, max_value = modulus1 - 1);
      uint32_t min BOOST_PREVENT_MACRO_SUBSTITUTION () const {return min_value;}
      uint32_t max BOOST_PREVENT_MACRO_SUBSTITUTION () const {return max_value;}

      static bool validation(uint32_t x) {
	return x == 831582319;
      }
    };

    template <typename Gen>
    class uniform_01_wrapper {
      Gen& gen;

      public:
      uniform_01_wrapper(Gen& gen): gen(gen) {}

      double operator()() const {return boost::uniform_01<double>()(gen);}

      typedef double result_type;
      BOOST_STATIC_CONSTANT(bool, has_fixed_range = true);
      static const double min_value, max_value;
      double min BOOST_PREVENT_MACRO_SUBSTITUTION () const {return min_value;}
      double max BOOST_PREVENT_MACRO_SUBSTITUTION () const {return max_value;}
    };

    template <typename Gen>
    const double uniform_01_wrapper<Gen>::min_value = 0.;
    template <typename Gen>
    const double uniform_01_wrapper<Gen>::max_value = 1.;
  }
}

#endif // BOOST_SPLITTABLE_ECUYER1988_HPP
