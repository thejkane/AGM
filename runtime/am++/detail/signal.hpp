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

#ifndef AMPLUSPLUS_DETAIL_SIGNAL_HPP
#define AMPLUSPLUS_DETAIL_SIGNAL_HPP

#include <vector>
#include <boost/smart_ptr.hpp>
#include <boost/assert.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/bind/apply.hpp>
#include <boost/ref.hpp>

namespace amplusplus {
  namespace detail {

template <typename HandlerBase, typename HandlerGen>
class signal_base {
  boost::shared_ptr<std::vector<boost::shared_ptr<HandlerBase> > > sig;

  public:
  signal_base(): sig(boost::make_shared<std::vector<boost::shared_ptr<HandlerBase> > >()) {}

  template <typename F>
  void* attach(const F& f) {
    boost::shared_ptr<HandlerBase> sp(new typename HandlerGen::template type<F>(f));
    boost::shared_ptr<std::vector<boost::shared_ptr<HandlerBase> > > new_p(boost::make_shared<std::vector<boost::shared_ptr<HandlerBase> > >(*sig));
    new_p->push_back(sp);
    sig.swap(new_p);
    return (void*)sp.get();
  }

  template <typename F>
  void* attach(const boost::reference_wrapper<F>& f) {
    boost::shared_ptr<HandlerBase> sp(new typename HandlerGen::template type<F&>(f));
    boost::shared_ptr<std::vector<boost::shared_ptr<HandlerBase> > > new_p(boost::make_shared<std::vector<boost::shared_ptr<HandlerBase> > >(*sig));
    new_p->push_back(sp);
    sig.swap(new_p);
    return (void*)sp.get();
  }

  void detach(void* handle) {
    boost::shared_ptr<std::vector<boost::shared_ptr<HandlerBase> > > new_p(boost::make_shared<std::vector<boost::shared_ptr<HandlerBase> > >(*sig));
    for (typename std::vector<boost::shared_ptr<HandlerBase> >::iterator
           i = new_p->begin(); i != new_p->end(); ++i) {
      if ((void*)(i->get()) == handle) {
        new_p->erase(i);
        sig.swap(new_p);
        return;
      }
    }
    BOOST_ASSERT (!"Detach of handle that was not found");
  }

  protected:
  template <typename Func>
  void call(const Func& func) const {
    boost::shared_ptr<std::vector<boost::shared_ptr<HandlerBase> > > sig_copy(sig);
    std::vector<boost::shared_ptr<HandlerBase> >& v(*sig_copy);
    for (typename std::vector<boost::shared_ptr<HandlerBase> >::iterator
           i = v.begin(); i != v.end(); ++i) {
      func(**i);
    }
  }
};

struct signal0_handler_base {
  virtual void operator()() = 0;
  virtual ~signal0_handler_base() {}
};

struct signal0_handler_gen {
  template <typename F>
  struct type: public signal0_handler_base {
    F f;
    type(typename boost::add_reference<typename boost::add_const<F>::type>::type f): f(f) {}
    void operator()() {f();}
  };
};

} // End namespace detail

class signal0: detail::signal_base<detail::signal0_handler_base, detail::signal0_handler_gen> {
  typedef detail::signal_base<detail::signal0_handler_base, detail::signal0_handler_gen> base_type;
  public:
  using base_type::attach;
  using base_type::detach;
  void operator()() {
    base_type::call(boost::bind(boost::apply<void>(), _1));
  }
};

namespace detail {

template <typename Arg>
struct signal1_handler_base {
  virtual void operator()(const Arg&) = 0;
  virtual ~signal1_handler_base() {}
};

template <typename Arg>
struct signal1_handler_gen {
  template <typename F>
  struct type: public signal1_handler_base<Arg> {
    F f;
    type(typename boost::add_reference<typename boost::add_const<F>::type>::type f): f(f) {}
    void operator()(const Arg& a) {f(a);}
  };
};

} // End namespace detail

template <typename Arg>
class signal1: detail::signal_base<detail::signal1_handler_base<Arg>, detail::signal1_handler_gen<Arg> > {
  typedef detail::signal_base<detail::signal1_handler_base<Arg>, detail::signal1_handler_gen<Arg> > base_type;
  public:
  using base_type::attach;
  using base_type::detach;
  void operator()(const Arg& a) {
    base_type::call(boost::bind(boost::apply<void>(), _1, boost::cref(a)));
  }
};

template <typename SigClass>
class scoped_attach: private boost::noncopyable {
  SigClass& sc;
  void* handle;

  public:
  template <typename F>
  scoped_attach(SigClass& sc, const F& f): sc(sc) {
    handle = sc.attach(f);
  }

  ~scoped_attach() {
    sc.detach(handle);
  }
};

}

#endif // AMPLUSPLUS_DETAIL_SIGNAL_HPP
