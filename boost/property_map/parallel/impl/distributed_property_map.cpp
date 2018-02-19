// Copyright (C) 2004-2006 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Nick Edmonds
//           Andrew Lumsdaine
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/bind.hpp>

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

namespace boost { namespace parallel {

PBGL_DISTRIB_PMAP_TEMPLATE
template<typename Reduce>
PBGL_DISTRIB_PMAP
::distributed_property_map(amplusplus::transport& trans, const GlobalMap& global,
                           const StorageMap& pm, const Reduce& reduce, 
			   MessageGenerator message_gen)
  : data(new data_t(trans, global, pm, reduce, reduce, 
		    Reduce::non_default_resolver, message_gen))
{
  data->ghost_cells.reset(new GhostCellStorage());
  data->reset = &data_t::template do_reset<Reduce>;
  data->ghost_cells->unset_limit();
  data->setup_messages();
}

PBGL_DISTRIB_PMAP_TEMPLATE
PBGL_DISTRIB_PMAP::~distributed_property_map() { }

PBGL_DISTRIB_PMAP_TEMPLATE
template<typename Reduce>
void 
PBGL_DISTRIB_PMAP::set_reduce(const Reduce& reduce)
{
  data->get_default_value = reduce;
  data->reduce = reduce;
  data->has_default_resolver = Reduce::non_default_resolver;
  int model = data->model;
  data->reset = &data_t::template do_reset<Reduce>;
  set_consistency_model(model);
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::prune_one_ghost_cell(const key_type& key, const value_type& value) const
{
  if (data->model & cm_flush) {
    // We need to flush values when we evict them.
    std::pair<key_type, value_type> const victim(key, value);
    data->put_msg.send(victim);
  }
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::prune_ghost_cells() const
{
  data->ghost_cells->prune(this);
}

PBGL_DISTRIB_PMAP_TEMPLATE
typename PBGL_DISTRIB_PMAP::value_type&
PBGL_DISTRIB_PMAP::cell(const key_type& key, bool request_if_missing) const
{
  bool found;
  value_type* value_ptr;
  data->ghost_cells->get_cell(key, value_ptr, found);

  value_type value;

  // Search for the ghost cell by key, and project back to the sequence
  if (!found) {
    if (data->has_default_resolver)
      // Since we have a default resolver, use it to create a default
      // value for this ghost cell.
      value = data->get_default_value(key);
    else { value = value_type(); }
      // TODO (NGE): Should we throw an exception when get() is called without a default resolver?

    typename PBGL_DISTRIB_PMAP::value_type& v = data->ghost_cells->insert_cell(key, value, this);

    if (request_if_missing)
      // Request the actual value of this key from its owner
      data->get_msg.send(std::make_pair(key,data->transport.rank()));
    return v;
  } else {
    return *value_ptr;
  }
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::data_t::
atomic_put_handler::operator() (const std::pair<key_type, value_type>& req) const
{
  using boost::get;

  owner_local_pair p = get(self->global, req.first);
  assert(p.first == self->transport.rank());

  put(self->storage, p.second, 
      self->reduce(req.first, get(self->storage, p.second), req.second));

//   value_type oldval, newval;

//   do {
//     oldval = get(self->storage, p.second);
//     newval = self->reduce(req.first, oldval, req.second);
//   } while(!maybe_exchange(self->storage, p.second, oldval, newval));
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::data_t::
put_handler::operator() (const std::pair<key_type, value_type>& req) const
{
  using boost::get;

  owner_local_pair p = get(self->global, req.first);
  assert(p.first == self->transport.rank());

  detail::maybe_put(self->storage, p.second,
                    self->reduce(req.first,
				 get(self->storage, p.second),
				 req.second));
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::data_t::
get_handler::operator() (const std::pair<key_type, rank_type>& msg) const
{
  using boost::get;

  rank_type source = msg.second;
  key_type key = msg.first;
  owner_local_pair p = get(self->global, key);
  value_type value = get(self->storage, p.second);

  self->get_reply_msg.send(boost::make_tuple(source, key, value));
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::data_t::
get_reply_handler::operator() (const boost::tuple<rank_type, key_type, value_type>& msg) const
{
  rank_type source = msg.template get<0>();
  key_type key = msg.template get<1>();
  value_type value = msg.template get<2>();

  value_type* value_ptr;
  bool found;
  self->ghost_cells->get_cell(key, value_ptr, found);
  if (found)
    *value_ptr = value;
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::data_t::
setup_messages()
{
  put_msg.set_handler(put_handler(*this));
  atomic_put_msg.set_handler(atomic_put_handler(*this));
  get_msg.set_handler(get_handler(*this));
  get_reply_msg.set_handler(get_reply_handler(*this));
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::set_consistency_model(int model)
{
  data->model = model;

  if (model & cm_backward) {
    // For backward consistency to work, we absolutely cannot throw
    // away any ghost cells.
    data->ghost_cells->unset_limit();
  }
}

PBGL_DISTRIB_PMAP_TEMPLATE
void
PBGL_DISTRIB_PMAP::set_max_ghost_cells(std::size_t max_ghost_cells)
{
  if ((data->model & cm_backward) && max_ghost_cells > 0)
      boost::throw_exception(std::runtime_error("distributed_property_map::set_max_ghost_cells: "
                                                "cannot limit ghost-cell usage with a backward "
                                                "consistency model"));

  data->ghost_cells->set_limit(max_ghost_cells, this);
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::clear()
{
  data->clear();
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::data_t::clear()
{
  ghost_cells->clear();
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::reset()
{
  if (data->reset) ((*data).*data->reset)();
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::flush()
{
  data->flush();
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::data_t::refresh_ghost_cells()
{
  using boost::get;

  // Collect the set of keys for which we will request values
  for (typename GhostCellStorage::iterator i = ghost_cells->begin();
       i != ghost_cells->end(); ++i)
    get_msg.send(std::make_pair(GhostCellStorage::get_key(i), transport.rank()));
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::data_t::flush()
{
  using boost::get;

  // Collect all of the flushed values
  for (typename GhostCellStorage::iterator i = ghost_cells->begin(); i != ghost_cells->end(); ++i)
    put_msg.send(std::make_pair(GhostCellStorage::get_key(i), GhostCellStorage::get_value(i)));
}

PBGL_DISTRIB_PMAP_TEMPLATE
void PBGL_DISTRIB_PMAP::synchronize()
{
  // Flush results
  if (data->model & cm_flush) {
    amplusplus::scoped_epoch epoch(data->transport);
    data->flush();
  }

  // Backward consistency
  if (data->model & cm_backward && !(data->model & (cm_clear | cm_reset))) {
    amplusplus::scoped_epoch epoch(data->transport);
    data->refresh_ghost_cells();
  }

  // Optionally clear results
  if (data->model & cm_clear)
    data->clear();

  // Optionally reset results
  if (data->model & cm_reset) {
    if (data->reset) ((*data).*data->reset)();
  }
}

PBGL_DISTRIB_PMAP_TEMPLATE
template<typename Resolver>
void PBGL_DISTRIB_PMAP::data_t::do_reset()
{
  Resolver* resolver = reduce.template target<Resolver>();
  assert(resolver);

  for (typename GhostCellStorage::iterator i = ghost_cells->begin(); i != ghost_cells->end(); ++i)
    GhostCellStorage::get_value(i) = (*resolver)(GhostCellStorage::get_key(i));
}

} } // end namespace boost::parallel
