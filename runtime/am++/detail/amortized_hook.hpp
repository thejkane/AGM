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

#ifndef AMPLUSPLUS_DETAIL_AMORTIZED_HOOK_HPP
#define AMPLUSPLUS_DETAIL_AMORTIZED_HOOK_HPP

#include <boost/intrusive/list.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread.hpp>
#include <vector>
#include <algorithm>
#include <am++/detail/append_buffer.hpp>
#include <am++/detail/thread_support.hpp>

#error "This code is not fixed up yet -- just use a Boost signal instead"

namespace amplusplus {
  namespace detail {

class amortized_hook_entry_base: public boost::intrusive::list_base_hook<> {
  public:
  virtual ~amortized_hook_entry_base() {}
  virtual void operator()() const = 0;
};

template <typename F, typename Arg>
class amortized_hook_entry: public amortized_hook_entry_base {
  public:
  amortized_hook_entry(amortized_hook& h, const F& f = F())
    : h(h), f(f)
  {
    h.entries.push_back(*this);
  }

  ~amortized_hook_entry() {
    h.entries.erase(h.entries.iterator_to(*this));
  }

  void* add(const Arg& a) { // FIXME: reuse empty slots
    slot_pair& loc = args.push_back(a);
    return (void*)(&loc);
  }

  void remove(void* loc_ptr) {
    ((slot_pair*)loc_ptr)->status = int(slot_empty);
  }

  virtual void operator()() const {
    for (typename append_buffer<slot_data>::const_iterator i = args.begin();
         i != args.end(); ++i) {
      if (i->status.load() == int(slot_full)) {
        f(i->data);
      }
    }
  }

  private:
  amortized_hook& h;
  F f;
  enum /* slot_status */ {slot_empty, slot_allocated, slot_full};
  struct slot_data {
    amplusplus::detail::atomic<int /* slot_status */> status;
    Arg data;
    slot_data(): status(int(slot_empty)), data() {}
    slot_data(const Arg& data): status(int(slot_full)), data(data) {}
    slot_data& operator=(const Arg& d) {
      this->status.store(int(slot_allocated)); // Mark that writing is in progress
      // Mem barrier here
      this->data = arg;
      // Mem barrier here
      this->status.store(int(slot_full)); // Now fully written
      return *this;
    }
  };
  append_buffer<slot_data> args;
};

class amortized_hook: boost::noncopyable {
  public:
  void operator()() const {
    for (boost::intrusive::list<amortized_hook_entry_base>::const_iterator
           i = entries.begin();
         i != entries.end(); ++i) {
      (*i)();
    }
  }
  private:
  boost::intrusive::list<amortized_hook_entry_base> entries;

  template <typename F, typename Arg>
  friend class amortized_hook_entry;
};

// Code modified from
// http://www.chaoticmind.net/~hcb/projects/boost.atomic/doc/atomic/usage_examples.html#boost_atomic.usage_examples.singleton
template <typename Tag>
class static_amortized_hook {
  public:
  static amortized_hook& get() {
    amortized_hook * tmp = instance_.load(boost::memory_order_consume);
    if (!tmp) {
      boost::mutex::scoped_lock l(instantiation_mutex);
      tmp=instance_.load(boost::memory_order_consume);
      if (!tmp) {
        tmp=new amortized_hook;
        instance_.store(tmp, boost::memory_order_release);
      }
    }
    return *tmp;
  }

  private:
  static amplusplus::detail::atomic<amortized_hook*> instance_;
  static boost::mutex instantiation_mutex;
};

  }
}

#endif // AMPLUSPLUS_DETAIL_AMORTIZED_HOOK_HPP
