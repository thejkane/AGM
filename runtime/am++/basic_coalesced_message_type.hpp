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

#ifndef AMPLUSPLUS_BASIC_COALESCED_MESSAGE_TYPE_HPP
#define AMPLUSPLUS_BASIC_COALESCED_MESSAGE_TYPE_HPP

#include <am++/traits.hpp>
#include <utility>
#include <boost/assert.hpp>
#include <vector>
#include <memory>
#include <boost/version.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/foreach.hpp>
#include <boost/noncopyable.hpp>
#include <am++/detail/typed_in_place_factory_owning.hpp>
#include <am++/performance_counters.hpp>
#include <stdint.h>
#include <am++/detail/signal.hpp>
#include <am++/detail/buffer_cache.hpp>
#include <am++/detail/thread_support.hpp>
#include <am++/detail/factory_wrapper.hpp>
#include <am++/dummy_buffer_sorter.hpp>

namespace amplusplus {

template <typename Arg, typename Handler, typename BufferSorter = amplusplus::DummyBufferSorter<Arg> >
class basic_coalesced_message_type;

struct basic_coalesced_message_type_gen {
  typedef transport::rank_type rank_type;
  size_t coalescing_size;

  public:
  template <typename Arg, typename Handler, typename BufferSorter = amplusplus::DummyBufferSorter<Arg> >
  struct inner {
    typedef basic_coalesced_message_type<Arg, Handler, BufferSorter> type;
  };

  explicit basic_coalesced_message_type_gen(size_t coalescing_size): coalescing_size(coalescing_size) {}
};

// Thread-safe functions:
//   Handler calls (from progress on underlying transport)
//   send
//   flush (both versions)
template <typename Arg, typename Handler, typename BufferSorter>
class basic_coalesced_message_type
#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
                                         : boost::noncopyable
#endif
{
  private:
  enum message_buffer_append_state {normal = 0, first_message = 1, buffer_now_full = 2, singleton_buffer = 3};

  struct message_buffer { // Not thread safe -- operations are locked by message type
    size_t count, max_count;
    boost::shared_ptr<void> data_owner;
    Arg* data;
    bool registered_with_td; // whether message_being_built has been called yet

    message_buffer(boost::shared_ptr<void> data_owner, size_t max_count)
        : count(0),
          max_count(max_count),
          data_owner(data_owner),
          data(static_cast<Arg*>(this->data_owner.get())),
          registered_with_td(false) {}

    message_buffer(): count(0), max_count(0), data_owner(), data(0), registered_with_td(false) {}

    void operator()() { // Deleter
      data = 0;
      data_owner.reset();
      registered_with_td = false;
    }

    bool valid() const {
      return data != 0;
    }

    message_buffer_append_state append(const Arg& x) {
      BOOST_ASSERT (count < max_count);
      int result = normal;
      if (count == 0) {
        result |= first_message;
      }
      data[count++] = x;
      if (count == max_count) result |= buffer_now_full;
      return message_buffer_append_state(result);
    }

    bool empty() const {
      BOOST_ASSERT (count < max_count);
      return count == 0;
    }

    void clear() {
      count = 0;
      registered_with_td = false;
    }

    size_t size() const {return count;}

    void swap(message_buffer& o) {
      std::swap(count, o.count);
      std::swap(max_count, o.max_count);
      data_owner.swap(o.data_owner);
      std::swap(data, o.data);
      std::swap(registered_with_td, o.registered_with_td);
    }
  };

  public:
  struct traits {
    typedef Arg arg_type;
    typedef Handler handler_type;
  };
  typedef Arg arg_type;
  typedef Handler handler_type;

  basic_coalesced_message_type(basic_coalesced_message_type_gen gen,
                               transport trans,
                               valid_rank_set possible_dests_ = valid_rank_set(),
                               valid_rank_set possible_sources_ = valid_rank_set(),
                               const BufferSorter& buf_sorter = DummyBufferSorter<Arg>())
    : trans(trans),
      mt(trans.create_message_type<Arg>()),
      buf_cache(new detail::buffer_cache(trans, gen.coalescing_size * sizeof(Arg))),
      handler(),
      outgoing_buffers(trans.size()),
      coalescing_size(gen.coalescing_size),
      lock(),
      buffer_sorter(buf_sorter)
  {
    if (!possible_dests_) possible_dests_ = boost::make_shared<detail::all_ranks>(trans.size());
    if (!possible_sources_) possible_sources_ = boost::make_shared<detail::all_ranks>(trans.size());
    this->mt.set_max_count(this->coalescing_size);
    this->mt.set_possible_dests(possible_dests_);
    this->mt.set_possible_sources(possible_sources_);
    typedef transport::rank_type rank_type;
    for (size_t i = 0; i < possible_dests_->count(); ++i) {
      rank_type r = possible_dests_->rank_from_index(i);
      BOOST_ASSERT (r < trans.size());
      outgoing_buffers[r] = alloc_buffer();
    }
  }

  ~basic_coalesced_message_type() {
    // These must be done in this order
    this->outgoing_buffers.clear();
    this->buf_cache.reset();
  }

#ifndef BOOST_NO_DEFAULTED_FUNCTIONS
  basic_coalesced_message_type(basic_coalesced_message_type&&) = default;
  basic_coalesced_message_type(const basic_coalesced_message_type&) = delete;
#endif

  private:
  void send_buffer(const message_buffer& buf, transport::rank_type dest) {
    // std::cerr << "send_buffer " << buf.count << " to " << dest << std::endl;
    assert (buf.registered_with_td);
    mt.send(buf.data, buf.count, dest, buf);
  }

  struct raw_message_handler {
    basic_coalesced_message_type& mt;
    raw_message_handler(basic_coalesced_message_type& _mt): mt(_mt) {}
    void operator()(transport::rank_type src, const Arg* buf, size_t count) const {
      amplusplus::performance_counters::hook_message_received(src, count, sizeof(Arg));
      handler_type& h = mt.get_handler();
      BufferSorter buf_sorter = mt.get_buffer_sorter();
      buf_sorter.sort(buf, buf + count);
      for (size_t i = 0; i < count; ++i) {
        h(src, buf[i]);
      }
    }
  };
  
  public:
  void set_handler(const handler_type& handler_) {
    handler.reset(new handler_type(handler_));
    this->mt.set_handler(raw_message_handler(*this));
  }

  BufferSorter& get_buffer_sorter() {
    return buffer_sorter;
  }
#ifndef BOOST_NO_RVALUE_REFERENCES
  void set_handler(handler_type&& handler_) {
    handler.reset(new handler_type(std::move(handler_)));
    this->mt.set_handler(raw_message_handler(*this));
  }
#endif

  handler_type& get_handler() {
    BOOST_ASSERT (handler);
    return *handler.get();
  }

  void send(const arg_type& arg, transport::rank_type dest) {
    BOOST_ASSERT (trans.is_valid_rank(dest));
    BOOST_ASSERT (this->mt.get_possible_dests()->is_valid(dest));
    boost::lock_guard<amplusplus::detail::recursive_mutex> my_lock(lock);
    // std::cerr << "send " << dest << std::endl;
    message_buffer& buf_ref = outgoing_buffers[dest];
    message_buffer_append_state s = buf_ref.append(arg);
    switch (s) {
      case normal: break;
      case first_message: {
        trans.add_flush_object(boost::bind(&basic_coalesced_message_type::flush, this, dest));
        if (!buf_ref.registered_with_td) {
          this->mt.message_being_built(dest);
          buf_ref.registered_with_td = true;
        }
        break;
      }
      case buffer_now_full:
      case singleton_buffer: {
        message_buffer buf = alloc_buffer();
        buf_ref.swap(buf);
        amplusplus::performance_counters::hook_full_buffer_send(dest, buf.size(), sizeof(Arg));
        if (!buf.registered_with_td) {
          this->mt.message_being_built(dest);
          buf.registered_with_td = true;
        }
        send_buffer(buf, dest);
        break;
      }
      default: abort();
    }
  }

  void send_with_tid(const arg_type& arg, transport::rank_type dest, int /*tid*/) {
    this->send(arg, dest);
  }

  void message_being_built(transport::rank_type dest) {
    BOOST_ASSERT (trans.is_valid_rank(dest));
    boost::lock_guard<amplusplus::detail::recursive_mutex> my_lock(lock);
    message_buffer& buf_ref = outgoing_buffers[dest];
    if (!buf_ref.registered_with_td) {
      buf_ref.registered_with_td = true;
      this->mt.message_being_built(dest);
    }
  }

  bool flush(transport::rank_type dest) {
    BOOST_ASSERT (trans.is_valid_rank(dest));
    // std::cerr << "flush " << dest << std::endl;
    boost::lock_guard<amplusplus::detail::recursive_mutex> my_lock(lock);
    message_buffer& buffer_to_send_ref = outgoing_buffers[dest];
    if (!buffer_to_send_ref.valid() || buffer_to_send_ref.empty()) return false;
    // std::cout << "flushing to " << dest << std::endl;
    message_buffer buffer_to_send = alloc_buffer();
    buffer_to_send.swap(buffer_to_send_ref);
    amplusplus::performance_counters::hook_flushed_message_size(dest, buffer_to_send.size(), sizeof(Arg));
    send_buffer(buffer_to_send, dest);
    return true;
  }

  transport get_transport() const {return trans;}

  private:
  transport trans;
  message_type<Arg> mt;
#if defined(BOOST_NO_CXX11_SMART_PTR) || BOOST_VERSION < 105000
  boost::scoped_ptr<detail::buffer_cache> buf_cache;
  boost::scoped_ptr<handler_type> handler;
#else
  std::unique_ptr<detail::buffer_cache> buf_cache;
  std::unique_ptr<handler_type> handler;
#endif
  std::vector<message_buffer> outgoing_buffers;
  size_t coalescing_size;
  amplusplus::detail::recursive_mutex lock;
  BufferSorter buffer_sorter;

  message_buffer alloc_buffer() {
    boost::shared_ptr<void> buf_owner = buf_cache->allocate();
    message_buffer buf(buf_owner, coalescing_size);
    return buf;
  }
};

}

#endif // AMPLUSPLUS_BASIC_COALESCED_MESSAGE_TYPE_HPP
