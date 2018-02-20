// Copyright 2017-2018 The Trustees of Indiana University.
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
//  Authors: Thejaka Kanelwala
//           Andrew Lumsdaine

#ifndef AMPLUSPLUS_SIZE_COALESCED_MESSAGE_TYPE_HPP
#define AMPLUSPLUS_SIZE_COALESCED_MESSAGE_TYPE_HPP

#include <am++/traits.hpp>
#include <boost/assert.hpp>
#include <utility>
#include <vector>
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

// This code needs to match promela/size_coalesced_message_type2.pr

typedef unsigned char COALESCE_TYPE;

namespace amplusplus {

template <typename Arg, typename Handler>
class size_coalesced_message_type;

struct size_coalesced_message_type_gen {
  typedef transport::rank_type rank_type;
  size_t coalescing_size;
  unsigned int priority;

  public:
  template <typename Arg, typename Handler>
  struct inner {
    typedef size_coalesced_message_type<Arg, Handler> type;
  };

  explicit size_coalesced_message_type_gen(size_t coalescing_size, int p=0): coalescing_size(coalescing_size), priority(p) {}
};



// Thread-safe functions:
//   Handler calls (from progress on underlying transport)
//   send
//   flush (both versions)
template <typename Arg, typename Handler>
class size_coalesced_message_type {
  private:
  struct sp_deleter {boost::shared_ptr<void> p; sp_deleter(const boost::shared_ptr<void>& p): p(p) {} void operator()() {p.reset();}};
  struct cbuf_deleter {COALESCE_TYPE* p; cbuf_deleter(COALESCE_TYPE*& pi): p(pi) {} void operator()() { delete [] p;}};

  struct message_buffer {
    amplusplus::detail::atomic<unsigned int> count_allocated, count_written;
    static const unsigned int sender_active = UINT_MAX - UINT_MAX / 2;
    static const unsigned int count_mask = sender_active - 1;
    unsigned int max_count;
    amplusplus::detail::atomic<bool> registered_with_td;
    boost::shared_ptr<void> data_owner;
    COALESCE_TYPE* data;
    // TODO : do i need to change this also ?
    struct {
      struct size_test {amplusplus::detail::atomic<unsigned int> a, b; unsigned int c; amplusplus::detail::atomic<bool> c2; boost::shared_ptr<void> d; Arg* e;};
      char padding[128 - sizeof(size_test)];
    } false_sharing_padding;

    explicit message_buffer(size_t max_count)
        : count_allocated(0),
          count_written(0),
          max_count(max_count),
          registered_with_td(false),
          data_owner(),
          data(0) {BOOST_ASSERT (max_count != 0);}

    message_buffer()
      : count_allocated(0), count_written(0), max_count(0), registered_with_td(false), data_owner(), data(0)
    {}

    message_buffer(const message_buffer& mb)
      : count_allocated(0),
        count_written(0),
        max_count(mb.max_count),
        registered_with_td(false),
        data_owner(mb.data_owner),
        data(mb.data)
    {
      // Only empty buffers should be copied
      BOOST_ASSERT (mb.count_allocated.load() == 0);
      BOOST_ASSERT (mb.count_written.load() == 0);
      BOOST_ASSERT (mb.registered_with_td.load() == false);
    }

    message_buffer& operator=(const message_buffer& mb) {
      // Only empty buffers should be copied, and only onto empty buffers
      BOOST_ASSERT (mb.count_allocated.load() == 0);
      BOOST_ASSERT (mb.count_written.load() == 0);
      BOOST_ASSERT (mb.registered_with_td.load() == false);
      BOOST_ASSERT (count_allocated.load() == 0);
      BOOST_ASSERT (count_written.load() == 0);
      BOOST_ASSERT (registered_with_td.load() == false);
      max_count = mb.max_count;
      data = mb.data;
      data_owner = mb.data_owner;
      return *this;
    }

    ~message_buffer() {
      BOOST_ASSERT (count_allocated.load() == 0);
      BOOST_ASSERT (count_written.load() == 0);
      BOOST_ASSERT (registered_with_td.load() == false);
    }

    bool empty() const {
      return count_allocated.load() == 0;
    }

    void clear(const boost::shared_ptr<void>& new_data_owner) {
      data_owner = new_data_owner;
      data = static_cast<COALESCE_TYPE*>(data_owner.get());
      registered_with_td.store(false);
      count_written.store(0u);
      count_allocated.store(0u); // This must be last since it releases other threads that may be waiting to add things to the buffer, and it must be at least a release to pick up the other writes
    }
  };

  public:
  struct traits {
    typedef Arg arg_type;
    typedef Handler handler_type;
  };
  typedef typename traits::arg_type arg_type;
  typedef typename traits::handler_type handler_type;

  size_coalesced_message_type(size_coalesced_message_type_gen gen,
                                 transport trans,
                                 valid_rank_set possible_dests_ = valid_rank_set(),
                                 valid_rank_set possible_sources_ = valid_rank_set())
    : trans(trans),
      mt(trans.create_message_type<COALESCE_TYPE>(gen.priority)),
      buf_cache(new detail::buffer_cache(trans, gen.coalescing_size * sizeof(COALESCE_TYPE))),
      handler(handler_type()),
      outgoing_buffers(trans.size()),
      last_active(trans.size()),
      coalescing_size(gen.coalescing_size),
      alive(boost::make_shared<bool>(true))
  {
    BOOST_ASSERT (coalescing_size != 0);
    if (!possible_dests_) possible_dests_ = boost::make_shared<detail::all_ranks>(trans.size());
    if (!possible_sources_) possible_sources_ = boost::make_shared<detail::all_ranks>(trans.size());
    this->mt.set_max_count(gen.coalescing_size);
    this->mt.set_handler(raw_message_handler(*this));
    this->mt.set_possible_dests(possible_dests_);
    this->mt.set_possible_sources(possible_sources_);
    typedef transport::rank_type rank_type;
    trans.add_flush_object(boost::bind(&size_coalesced_message_type::flush, this, alive));
    for (size_t i = 0; i < possible_dests_->count(); ++i) {
      rank_type r = possible_dests_->rank_from_index(i);
      BOOST_ASSERT (r < trans.size());
      outgoing_buffers[r] = message_buffer(gen.coalescing_size);
      outgoing_buffers[r].clear(buf_cache->allocate());
      last_active[r] = 0;
    }
  }

#ifndef BOOST_NO_DEFAULTED_FUNCTIONS
  size_coalesced_message_type(size_coalesced_message_type&& o)
    : trans(o.trans), mt(std::move(o.mt)), buf_cache(),
      handler(std::move(o.handler)), outgoing_buffers(std::move(o.outgoing_buffers)),
      last_active(std::move(o.last_active)),
      coalescing_size(o.coalescing_size) {
    buf_cache.swap(o.buf_cache);
  }
#endif

  ~size_coalesced_message_type() {
    // These must be done in this order
    *alive = false;
    this->outgoing_buffers.clear();
    this->buf_cache.reset();
  }

  private:
  bool send_buffer(message_buffer& buf, unsigned int my_id, transport::rank_type dest) {
    // fprintf(stderr, "%zu: send_buffer %08x to %zu, current state is %08x\n", trans.rank(), my_id, dest, buf.count_allocated.load());
    BOOST_ASSERT (buf.count_allocated.load() & message_buffer::sender_active);
    const unsigned int count = my_id & message_buffer::count_mask;
    if ((my_id & message_buffer::sender_active) != 0) return false;
    BOOST_ASSERT (!!buf.data_owner);
    BOOST_ASSERT (count <= buf.max_count);
    if (count != 0) {
      while (buf.count_written.load() != count) {amplusplus::detail::do_pause();} // Make sure all elements are written, must be an acquire so we can read data in the buffer
      BOOST_ASSERT (buf.registered_with_td.load() == true); // This may not be true before the previous line's while loop is done
      COALESCE_TYPE* send_data = buf.data;
      boost::shared_ptr<void> send_data_owner = buf.data_owner;
      buf.clear(buf_cache->allocate());
      // fprintf(stderr, "%zu: actual send %u to %zu, current state is %08x\n", trans.rank(), count, dest, buf.count_allocated.load());
      this->mt.send(send_data, count, dest, sp_deleter(send_data_owner));
      return true;
    }
    return false;
  }
 
  
  struct raw_message_handler {
    size_coalesced_message_type& mt;
    raw_message_handler(size_coalesced_message_type& mt): mt(mt) {}
    void operator()(transport::rank_type src, COALESCE_TYPE* buf, size_t count) const {
      amplusplus::performance_counters::hook_message_received(src, count, sizeof(COALESCE_TYPE));
      handler_type& h = mt.handler;
      size_t currentpos = 0;
      while(count != 0) {
	arg_type at;
	size_t sz;
	std::memcpy(&sz, buf+currentpos, szsize);
	assert(sz > szsize);
	at.deserialize(buf+currentpos+szsize, (sz-szsize));
	count -= sz; // sz includes szsize
	currentpos += sz;
	h(src, at);
      }
    }
  };
  

  public:
  void set_handler(const handler_type& handler_) {
    handler = handler_;
  }

#ifndef BOOST_NO_RVALUE_REFERENCES
  void set_handler(handler_type&& handler_) {handler = std::move(handler_);}
#endif

  handler_type& get_handler() {
    return handler;
  }


  void send(arg_type& wi, transport::rank_type dest) {
    BOOST_ASSERT (trans.is_valid_rank(dest));

    size_t sz = wi.get_size() + szsize;
    message_buffer& buf = outgoing_buffers[dest];
    // if current message is greater than the max_count send the message
    // need to be cautios on this send -- are we deleting the databuf ?
    if (sz > buf.max_count) {
      fprintf(stderr, "[ERROR] Sending message size is greater than the the coalescing size. Please increase the coalescing size at least upto %zu.\n", sz);
      abort();
      return;
    }

    const unsigned int max_count = buf.max_count;
    while (true) {
      unsigned int oldval = 0;
      while (true) {
        assert (max_count > 0);
        oldval = buf.count_allocated.load();
        // fprintf(stderr, "send_spin %08x of %x\n", x, max_count);
        if ((oldval & message_buffer::count_mask) < max_count && (oldval & message_buffer::sender_active) == 0) {
          break;
        }
        amplusplus::detail::do_pause();
      }

      unsigned int my_id = (oldval+sz);
      if (buf.count_allocated.compare_exchange_strong(oldval, my_id)) {
	// is this a new buffer, then register with td
	if (oldval == 0) {
	  if (!buf.registered_with_td.exchange(true)) {
	    //	    fprintf(stderr, ":register:");
	    this->mt.message_being_built(dest);
	  }
	}

	// only one thread should come here
	if ((my_id & message_buffer::count_mask) >= max_count) { 
	  buf.count_allocated.store(message_buffer::sender_active);
	  // send the buffer
	  amplusplus::performance_counters::hook_full_buffer_send(dest, oldval, sizeof(COALESCE_TYPE));
	  send_buffer(buf, oldval, dest);

	  // continue again
	  continue;
	}
      } else
	continue;

      BOOST_ASSERT ((my_id & message_buffer::count_mask) < max_count);
      BOOST_ASSERT (buf.data_owner);
      BOOST_ASSERT (buf.data);
      //      fprintf(stderr, ",send2");
      // put size
      std::memcpy(buf.data+oldval, &sz, szsize);
      wi.serialize(buf.data+oldval+szsize);
      buf.count_written.fetch_add(sz);
      break;
    }
  }

  void send_with_tid(const arg_type& arg, transport::rank_type dest, int /*tid*/) {
    this->send(arg, dest);
  }

  void message_being_built(transport::rank_type dest) {
    BOOST_ASSERT (trans.is_valid_rank(dest));
    message_buffer& buf = this->outgoing_buffers[dest];
    bool was_registered = buf.registered_with_td.exchange(true);
    if (!was_registered) {
      this->mt.message_being_built(dest);
    }
  }

  bool flush(boost::shared_ptr<bool> alive) {
    //std::cout << "inside flush" << std::endl;
    BOOST_ASSERT (alive.get());
    if(!*alive) return false;
    //    fprintf(stderr, "%zu: size_coalesced_message_type::flush\n", trans.rank());
    for (size_t i = 0; i < this->mt.get_possible_dests()->count(); ++i) {
      rank_type r = this->mt.get_possible_dests()->rank_from_index(i);
      BOOST_ASSERT (trans.is_valid_rank(r));
      message_buffer& buf = this->outgoing_buffers[r];
      const unsigned int max_count = buf.max_count;
      unsigned int my_id = buf.count_allocated.load();
      if(my_id != last_active[r].load()) {
	last_active[r].store(my_id);
	continue;
      }
      while (true) {
	if (my_id > 0 && my_id < max_count) {
	  if (buf.count_allocated.compare_exchange_weak(my_id, message_buffer::sender_active)) break;
	  amplusplus::detail::do_pause();
	} else {
	  break;
	}
      }
      if (my_id > 0 && my_id < max_count) {
	amplusplus::performance_counters::hook_flushed_message_size(r, my_id, sizeof(COALESCE_TYPE));
	//	fprintf(stderr, "sending the buffer");
	send_buffer(buf, my_id, r);
      }
    }
    return true;
  }

  transport get_transport() const {return trans;}

  private:
  transport trans;
  message_type<COALESCE_TYPE> mt;
  boost::scoped_ptr<detail::buffer_cache> buf_cache;
  handler_type handler;
  std::vector<message_buffer> outgoing_buffers;
  std::vector<amplusplus::detail::atomic<unsigned int> > last_active;
  size_t coalescing_size;
  boost::shared_ptr<bool> alive;
  static const size_t szsize = sizeof(uint64_t);
};

}

#endif // AMPLUSPLUS_SIZE_COALESCED_MESSAGE_TYPE_HPP
