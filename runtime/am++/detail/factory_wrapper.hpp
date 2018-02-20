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

#ifndef AMPLUSPLUS_DETAIL_FACTORY_WRAPPER_HPP
#define AMPLUSPLUS_DETAIL_FACTORY_WRAPPER_HPP

#include <boost/preprocessor.hpp>
#include <am++/detail/typed_in_place_factory_owning.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/utility/enable_if.hpp>

#define AMPLUSPLUS_MAKE_WRAPPER(type, tparam_seq) \
AMPLUSPLUS_ADD_WRAP_IF_NONZERO(BOOST_PP_SEQ_SIZE(tparam_seq))(template <, >, BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AMPLUSPLUS_PREPEND_TYPENAME, no_data, tparam_seq))) \
class type { \
  public: \
  typedef BOOST_PP_CAT(type, _cls)AMPLUSPLUS_ADD_WRAP_IF_NONZERO(BOOST_PP_SEQ_SIZE(tparam_seq))(<, >, BOOST_PP_SEQ_ENUM(tparam_seq)) underlying_type; \
  typedef amplusplus::message_type_traits<underlying_type> traits; \
  AMPLUSPLUS_MOVE_CONSTRUCTOR_MAYBE(type) \
  AMPLUSPLUS_MAKE_WRAPPER_BODY(type) \
};

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
#define AMPLUSPLUS_MOVE_CONSTRUCTOR_MAYBE(type) /**/
#else
#define AMPLUSPLUS_MOVE_CONSTRUCTOR_MAYBE(type) type(type&&) = default;
#endif

#define AMPLUSPLUS_ADD_WRAP_IF_NONZERO(n) \
  BOOST_PP_IF(n, AMPLUSPLUS_ACTUALLY_ADD_WRAP, BOOST_PP_TUPLE_EAT(3))
#define AMPLUSPLUS_ACTUALLY_ADD_WRAP(left, right, body) left body right

#define AMPLUSPLUS_PREPEND_TYPENAME(s, ignored, tn) typename tn

#define AMPLUSPLUS_MAKE_WRAPPER_BODY(type_) \
  public: \
  underlying_type underlying; \
 \
  type_(): underlying() {} \
  template <typename A0> type_(const A0& a0, typename boost::disable_if<boost::is_base_and_derived<amplusplus::detail::typed_in_place_factory_owning_base, A0>, void*>::type = 0): underlying(a0) {} \
  BOOST_PP_REPEAT_FROM_TO(2, 11, AMPLUSPLUS_ONE_FORWARDING_CONSTRUCTOR, type_) \
  type_(const amplusplus::detail::typed_in_place_factory_owning0<underlying_type>&): underlying() {} \
  template <typename A0> type_(const amplusplus::detail::typed_in_place_factory_owning1<underlying_type, A0>& f): underlying(f.a0) {} \
  BOOST_PP_REPEAT_FROM_TO(2, 11, AMPLUSPLUS_ONE_FACTORY_CONSTRUCTOR, type_) \
 \
  underlying_type& get() {return underlying;} \
  const underlying_type& get() const {return underlying;} \
  underlying_type* operator->() {return &underlying;} \
  const underlying_type* operator->() const {return &underlying;}

#define AMPLUSPLUS_ONE_FORWARDING_CONSTRUCTOR(z, i, type) \
  template <BOOST_PP_ENUM_PARAMS(i, typename A)> type(BOOST_PP_ENUM_BINARY_PARAMS(i, const A, & a)): underlying(BOOST_PP_ENUM_PARAMS(i, a)) {}

#define AMPLUSPLUS_ONE_FACTORY_CONSTRUCTOR(z, i, type) \
  template <BOOST_PP_ENUM_PARAMS(i, typename A)> type(const amplusplus::detail::BOOST_PP_CAT(typed_in_place_factory_owning, i)<underlying_type BOOST_PP_ENUM_TRAILING_PARAMS(i, A)>& f): underlying(BOOST_PP_ENUM_PARAMS(i, f.a)) {}

#endif // AMPLUSPLUS_DETAIL_FACTORY_WRAPPER_HPP
