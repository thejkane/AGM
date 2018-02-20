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

#ifndef AMPLUSPLUS_DETAIL_TYPE_INFO_MAP_HPP
#define AMPLUSPLUS_DETAIL_TYPE_INFO_MAP_HPP

#include <typeinfo>
#include <map>
#include <utility>

namespace amplusplus {
  namespace detail {

template <typename T>
struct type_info_wrapper {virtual ~type_info_wrapper() {}}; // To enable RTTI

template <typename T>
class type_info_map {
  public:
  type_info_map(): m() {}

  void insert(const std::type_info& ti, const T& val) {
    m.insert(std::make_pair(&ti, val));
  }

  const T* lookup(const std::type_info& ti) {
    typename underlying_type::const_iterator i = m.find(&ti);
    if (i == m.end()) return 0;
    return &i->second;
  }

  void clear() {m.clear();}

  private:
  struct type_info_compare {
    bool operator()(const std::type_info* a, const std::type_info* b) const {
      return a->before(*b);
    }
  };
  typedef std::map<const std::type_info*, T, type_info_compare> underlying_type;
  underlying_type m;
};

template <typename T>
const std::type_info& get_type_info() {return typeid(type_info_wrapper<T>);}

  }
}

#endif // AMPLUSPLUS_DETAIL_TYPE_INFO_MAP_HPP
