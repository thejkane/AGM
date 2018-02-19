// Copyright 2005-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine
#ifndef BOOST_PARALLEL_GLOBAL_INDEX_MAP_HPP
#define BOOST_PARALLEL_GLOBAL_INDEX_MAP_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include "am++/basic_coalesced_message_type.hpp"

#include <boost/property_map/property_map.hpp>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace boost { namespace parallel {

template<typename IndexMap, typename GlobalMap>
class global_index_map
{
public:
  typedef typename property_traits<IndexMap>::key_type key_type;
  typedef typename property_traits<IndexMap>::value_type value_type;
  typedef value_type reference;
  typedef readable_property_map_tag category;

private:
  typedef typename amplusplus::transport::rank_type rank_type;

private:

  // Message handlers
  struct send_count_handler {
    
    send_count_handler() : self(NULL) {}
    send_count_handler(global_index_map* self) : self(self) {}

    void operator() (int source, const value_type& num_local_indices)
    {
      self->starting_index[source + 1] = num_local_indices;
    }

  protected:
    global_index_map* self;
  };
  friend struct send_count_handler;

  struct set_starting_index_handler {

    set_starting_index_handler() : self(NULL) {}
    set_starting_index_handler(global_index_map* self) : self(self) {}

    void operator() (int /* source */, const std::pair<rank_type, value_type>& data)
    {
      self->starting_index[data.first] = data.second;
    }

  protected:
    global_index_map* self;
  };
  friend struct set_starting_index_handler;

  // Message types
  typedef amplusplus::basic_coalesced_message_type<value_type, send_count_handler> 
      send_count_message_type;

  typedef amplusplus::basic_coalesced_message_type<std::pair<rank_type, value_type>,
						   set_starting_index_handler> 
      set_starting_index_message_type;

public:
  global_index_map(amplusplus::transport& trans, value_type num_local_indices, 
                   IndexMap index_map, GlobalMap global)
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<value_type>(),
      amplusplus::register_mpi_datatype<std::pair<rank_type, value_type> >(), 0)),
      starting_index(trans.size() + 1), index_map(index_map), global(global),
      send_count_msg(amplusplus::basic_coalesced_message_type_gen(1 << 12), trans), 
      set_starting_index_msg(amplusplus::basic_coalesced_message_type_gen(1 << 12), trans)
  {
    send_count_msg.set_handler(send_count_handler(this));
    set_starting_index_msg.set_handler(set_starting_index_handler(this));

    {
      amplusplus::scoped_epoch epoch(trans);

      send_count_msg.send(num_local_indices, 0);
    }
    
    starting_index[0] = 0;
    {
      amplusplus::scoped_epoch epoch(trans);

      if (trans.rank() == 0) {
	
	// Calculate starting indices from received counts
	for (rank_type src = 0; src < trans.size(); ++src)
	  starting_index[src + 1] += starting_index[src];

	// Broadcast calculated indices vector
	// TODO (NGE): If we had serialization/datatypes for vectors we could eliminate the inner loop
	for (rank_type dest = 1; dest < trans.size(); ++dest)
	  for (rank_type idx = 1; idx < trans.size() + 1; ++idx)
	    set_starting_index_msg.send(std::make_pair(idx, starting_index[idx]), dest);
      }
    }
  }

  friend inline value_type 
  get(const global_index_map& gim, const key_type& x)
  {
    using boost::get;
    return gim.starting_index[get(gim.global, x).first]
           + get(gim.index_map, x);
  }

private:
  const int dummy_first_member_for_init_order; // Unused

  std::vector<value_type> starting_index;
  IndexMap index_map;
  GlobalMap global;

  send_count_message_type send_count_msg;
  set_starting_index_message_type set_starting_index_msg;
};

} } // end namespace boost::parallel

#endif // BOOST_PARALLEL_GLOBAL_INDEX_MAP_HPP
