// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#include <boost/optional.hpp>
#include <cassert>
#include <functional>
#include <algorithm>

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

namespace boost { namespace graph { namespace distributed {

template<BOOST_DISTRIBUTED_QUEUE_PARMS>
BOOST_DISTRIBUTED_QUEUE_TYPE::
distributed_queue(amplusplus::transport& trans, const OwnerMap& owner,
                  boost::shared_ptr<Buffer> buffer, const UnaryPredicate& pred, 
		  MessageGenerator message_gen)
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<value_type>(), 0)),
    transport(trans),
    owner(owner),
    buffer(buffer),
    pred(pred),
    push_msg(message_gen, transport, owner, amplusplus::duplicate_removal(boost::parallel::identity<value_type>()))
{
  push_msg.set_handler(push_handler(*this));
}

template<BOOST_DISTRIBUTED_QUEUE_PARMS>
BOOST_DISTRIBUTED_QUEUE_TYPE::
distributed_queue(amplusplus::transport& trans, const OwnerMap& owner,
                  boost::shared_ptr<Buffer> buffer, boost::shared_ptr<Buffer> incoming_buffer, 
		  const UnaryPredicate& pred, MessageGenerator message_gen)
  : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<value_type>(), 0)),
    transport(trans),
    owner(owner),
    buffer(buffer),
    incoming_buffer(incoming_buffer),
    pred(pred),
    push_msg(message_gen, trans, owner, amplusplus::duplicate_removal(boost::parallel::identity<value_type>()))
{
  push_msg.set_handler(push_handler(*this));
}

template<BOOST_DISTRIBUTED_QUEUE_PARMS>
struct BOOST_DISTRIBUTED_QUEUE_TYPE::push_handler {
  
  push_handler() : self(NULL) {}
  push_handler(distributed_queue& self) : self(&self) {}
  
  void operator() (const value_type& data) const
  {
    if (self->pred(data)) {
      if (self->incoming_buffer)
	self->incoming_buffer->push_back(data);
      else
	self->buffer->push_back(data);
    }
  }
  
protected:
  distributed_queue* self;
};

} } } // end namespace boost::graph::distributed
