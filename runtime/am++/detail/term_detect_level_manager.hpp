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

#ifndef AMPLUSPLUS_DETAIL_TERM_DETECT_LEVEL_MANAGER_HPP
#define AMPLUSPLUS_DETAIL_TERM_DETECT_LEVEL_MANAGER_HPP

#include <set>
#include <functional>
#include <boost/thread.hpp>
#include <boost/assert.hpp>
#include <am++/detail/thread_support.hpp>

namespace amplusplus {
  namespace detail {

class term_detect_level_manager {
  std::multiset<size_t, std::greater<size_t> > td_levels_requested;
  size_t current_td_level;
  mutable amplusplus::detail::recursive_mutex my_lock;

  void update_td_level() { // my_lock must be held
    current_td_level = td_levels_requested.empty() ?
                       0 :
                       *td_levels_requested.begin();
  }

  public:
  size_t get() const {return current_td_level;}

  void insert(size_t l) {
    boost::lock_guard<amplusplus::detail::recursive_mutex> lock(my_lock);
    td_levels_requested.insert(l);
    update_td_level();
  }

  void erase(size_t l) {
    boost::lock_guard<amplusplus::detail::recursive_mutex> lock(my_lock);
    BOOST_ASSERT (td_levels_requested.find(l) != td_levels_requested.end());
    td_levels_requested.erase(td_levels_requested.find(l));
    update_td_level();
  }
};

  }
}

#endif // AMPLUSPLUS_DETAIL_TERM_DETECT_LEVEL_MANAGER_HPP
