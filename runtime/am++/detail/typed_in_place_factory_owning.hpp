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

#ifndef AMPLUSPLUS_DETAIL_TYPED_IN_PLACE_FACTORY_OWNING_HPP
#define AMPLUSPLUS_DETAIL_TYPED_IN_PLACE_FACTORY_OWNING_HPP

#include <boost/preprocessor.hpp>
#include <boost/utility/typed_in_place_factory.hpp>

// Like boost::typed_in_place_factoryN, but owning the data used for
// construction to avoid dangling references.

namespace amplusplus {
  namespace detail {
    // This inheritance is to trick Boost.Optional into allowing our factory as a valid argument type
    struct typed_in_place_factory_owning_base: boost::typed_in_place_factory_base {};

#define AMPLUSPLUS_DETAIL_MAKE_ONE_TIPFO(repdim, nparams, data) \
    template <typename Obj BOOST_PP_ENUM_TRAILING_PARAMS_Z(repdim, nparams, typename A)> \
    class BOOST_PP_CAT(typed_in_place_factory_owning, nparams): public typed_in_place_factory_owning_base { \
      typedef boost::BOOST_PP_CAT(typed_in_place_factory, nparams)<Obj BOOST_PP_ENUM_TRAILING_PARAMS_Z(repdim, nparams, A)> base_type; \
      \
      public: \
      typedef Obj value_type; \
      explicit BOOST_PP_CAT(typed_in_place_factory_owning, nparams)(BOOST_PP_ENUM_BINARY_PARAMS_Z(repdim, nparams, const A, & a)) BOOST_PP_IF(nparams, :, /**/) BOOST_PP_CAT(BOOST_PP_ENUM_, repdim)(nparams, AMPLUSPLUS_DETAIL_MAKE_ONE_TIPFO_INIT, x) {} \
      void apply(void* obj) const {base_type(BOOST_PP_ENUM_PARAMS_Z(repdim, nparams, a)).apply(obj);} \
      \
      public: \
      BOOST_PP_CAT(BOOST_PP_REPEAT_, repdim)(nparams, AMPLUSPLUS_DETAIL_MAKE_ONE_TIPFO_MEM, x) \
    };

#define AMPLUSPLUS_DETAIL_MAKE_ONE_TIPFO_INIT(repdim, i, data) BOOST_PP_CAT(a, i)(BOOST_PP_CAT(a, i))
#define AMPLUSPLUS_DETAIL_MAKE_ONE_TIPFO_MEM(repdim, i, data) BOOST_PP_CAT(A, i) BOOST_PP_CAT(a, i);

    BOOST_PP_REPEAT(11, AMPLUSPLUS_DETAIL_MAKE_ONE_TIPFO, x)
  }
}

#endif
