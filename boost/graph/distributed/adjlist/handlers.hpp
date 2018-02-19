// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds

// This file contains code for the distributed adjacency list's
// message handlers. It should not be included directly by users.

#ifndef BOOST_GRAPH_DISTRIBUTED_ADJLIST_HANDLERS_HPP
#define BOOST_GRAPH_DISTRIBUTED_ADJLIST_HANDLERS_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/graph/parallel/detail/untracked_pair.hpp>

namespace boost {

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
void
PBGL_DISTRIB_ADJLIST_TYPE::
setup_messages()
{
  add_edge_msg.set_handler(add_edge_handler(*this));
  add_edge_with_property_msg.set_handler(add_edge_with_property_handler(*this));
  remove_edge_msg.set_handler(remove_edge_handler(*this));
}
    
template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
struct PBGL_DISTRIB_ADJLIST_TYPE::add_edge_handler {

  add_edge_handler() : self(NULL) {}
  add_edge_handler(adjacency_list& self) : self(&self) {}

  void operator() (int /* source */, const msg_add_edge_data& data)
  {
    add_edge(vertex_descriptor(self->processor(), data.first), 
	     data.second, *self);

    std::cout << self->processor() << ": remote add edge " << data.first << "@" 
	      << self->processor() << " -> " << data.second.local << "@" << data.second.owner
	      << std::endl;
  }

protected:
  adjacency_list* self;
};

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
struct PBGL_DISTRIB_ADJLIST_TYPE::add_edge_with_property_handler {

  add_edge_with_property_handler() : self(NULL) {}
  add_edge_with_property_handler(adjacency_list& self) : self(&self) {}

  void operator() (int /* source */, const msg_add_edge_with_property_data& data)
  {
    add_edge(vertex_descriptor(self->processor(), get<0>(data)), 
	     get<1>(data), get<2>(data), *self);
  }

protected:
  adjacency_list* self;
};

template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
struct PBGL_DISTRIB_ADJLIST_TYPE::remove_edge_handler {

  remove_edge_handler() : self(NULL) {}
  remove_edge_handler(adjacency_list& self) : self(&self) {}

  void operator() (int /* source */, const msg_remove_edge_data& )
  {
    // TODO (NGE): Define appropriate operation here
  }

protected:
  adjacency_list* self;
};

} // end namespace boost

#endif // BOOST_GRAPH_DISTRIBUTED_ADJLIST_HANDLERS_HPP

