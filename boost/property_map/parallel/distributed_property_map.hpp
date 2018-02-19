// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

// The placement of this #include probably looks very odd relative to
// the #ifndef/#define pair below. However, this placement is
// extremely important to allow the various property map headers to be
// included in any order.
#include <boost/property_map/property_map.hpp>

#ifndef BOOST_PARALLEL_DISTRIBUTED_PROPERTY_MAP_HPP
#define BOOST_PARALLEL_DISTRIBUTED_PROPERTY_MAP_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

// #define THREADED_DPMAP

#include "am++/counter_coalesced_message_type.hpp"

#include <boost/property_map/property_map.hpp>

#include <boost/graph/parallel/graph_utility.hpp> // for map_of_project1stdot1st
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/graph/detail/edge.hpp>
#include <boost/function/function1.hpp>
#include <boost/function/function3.hpp>
#include <vector>
#include <set>
#include <boost/graph/parallel/basic_reduce.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/property_map/parallel/local_property_map.hpp>
#include <map>
#include <boost/version.hpp>
#include <boost/property_map/parallel/impl/multi_index_ghost_cell.hpp>
#include <boost/property_map/parallel/impl/associative_cache_ghost_cell.hpp>
#include <boost/property_map/parallel/impl/unordered_map_ghost_cell.hpp>
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap

namespace boost { namespace parallel {

namespace detail {
  /**************************************************************************
   * Metafunction that degrades an Lvalue Property Map category tag to
   * a Read Write Property Map category tag.
   **************************************************************************/
  template<bool IsLvaluePropertyMap>
  struct make_nonlvalue_property_map
  {
    template<typename T> struct apply { typedef T type; };
  };

  template<>
  struct make_nonlvalue_property_map<true>
  {
    template<typename>
    struct apply
    {
      typedef read_write_property_map_tag type;
    };
  };

  /**************************************************************************
   * Performs a "put" on a property map so long as the property map is
   * a Writable Property Map or a mutable Lvalue Property Map. This
   * is required because the distributed property map's message
   * handler handles "put" messages even for a const property map,
   * although receipt of a "put" message is ill-formed.
   **************************************************************************/
  template<bool IsLvaluePropertyMap>
  struct maybe_put_in_lvalue_pm
  {
    template<typename PropertyMap, typename Key, typename Value>
    static inline void
    do_put(PropertyMap, const Key&, const Value&)
    { assert(false); }
  };

  template<>
  struct maybe_put_in_lvalue_pm<true>
  {
    template<typename PropertyMap, typename Key, typename Value>
    static inline void
    do_put(PropertyMap pm, const Key& key, const Value& value)
    { 
      using boost::put;

      put(pm, key, value); 
    }
  };

  template<typename PropertyMap, typename Key, typename Value>
  inline void
  maybe_put_impl(PropertyMap pm, const Key& key, const Value& value,
                 writable_property_map_tag)
  {
    using boost::put;

    put(pm, key, value);
  }

  template<typename PropertyMap, typename Key, typename Value>
  inline void
  maybe_put_impl(PropertyMap pm, const Key& key, const Value& value,
                 lvalue_property_map_tag)
  {
    typedef typename property_traits<PropertyMap>::value_type value_type;
    typedef typename property_traits<PropertyMap>::reference reference;
    // DPG TBD: Some property maps are improperly characterized as
    // lvalue_property_maps, when in fact they do not provide true
    // references. The most typical example is those property maps
    // built from vector<bool> and its iterators, which deal with
    // proxies. We don't want to mischaracterize these as not having a
    // "put" operation, so we only consider an lvalue_property_map as
    // constant if its reference is const value_type&. In fact, this
    // isn't even quite correct (think of a
    // vector<bool>::const_iterator), but at present C++ doesn't
    // provide us with any alternatives.
    typedef is_same<const value_type&, reference> is_constant;

    maybe_put_in_lvalue_pm<(!is_constant::value)>::do_put(pm, key, value);
  }

  template<typename PropertyMap, typename Key, typename Value>
  inline void
  maybe_put_impl(PropertyMap, const Key&, const Value&, ...)
  { assert(false); }

  template<typename PropertyMap, typename Key, typename Value>
  inline void
  maybe_put(PropertyMap pm, const Key& key, const Value& value)
  {
    maybe_put_impl(pm, key, value,
                   typename property_traits<PropertyMap>::category());
  }
} // end namespace detail

/** The consistency model used by the distributed property map. */
enum consistency_model {
  cm_forward = 1 << 0,
  cm_backward = 1 << 1,
  cm_bidirectional = cm_forward | cm_backward,
  cm_flush = 1 << 2,
  cm_reset = 1 << 3,
  cm_clear = 1 << 4
};

/** Distributed property map adaptor.
 *
 *  The distributed property map adaptor is a property map whose
 *  stored values are distributed across multiple non-overlapping
 *  memory spaces on different processes. Values local to the current
 *  process are stored within a local property map and may be
 *  immediately accessed via @c get and @c put. Values stored on
 *  remote processes may also be access via @c get and @c put, but the
 *  behavior differs slightly:
 *
 *  - @c put operations update a local ghost cell and send a "put"
 *    message to the process that owns the value. The owner is free to
 *    update its own "official" value or may ignore the put request.
 *
 *  - @c get operations returns the contents of the local ghost
 *    cell. If no ghost cell is available, one is created using the
 *    default value provided by the "reduce" operation. See, e.g.,
 *    @ref basic_reduce and @ref property_reduce.
 *
 * Using distributed property maps requires a bit more care than using
 * local, sequential property maps. While the syntax and semantics are
 * similar, distributed property maps may contain out-of-date
 * information that can only be guaranteed to be synchronized by
 * calling the @ref synchronize function in all processes.
 *
 * To address the issue of out-of-date values, distributed property
 * maps are supplied with a reduction operation. The reduction
 * operation has two roles:
 *
 *   -# When a value is needed for a remote key but no value is
 *      immediately available, the reduction operation provides a
 *      suitable default. For instance, a distributed property map
 *      storing distances may have a reduction operation that returns
 *      an infinite value as the default, whereas a distributed
 *      property map for vertex colors may return white as the
 *      default.
 *
 *   -# When a value is received from a remote process, the process
 *      owning the key associated with that value must determine which
 *      value---the locally stored value, the value received from a
 *      remote process, or some combination of the two---will be
 *      stored as the "official" value in the property map. The
 *      reduction operation transforms the local and remote values
 *      into the "official" value to be stored.
 *
 * @tparam StorageMap the type of the property map that will
 * store values for keys local to this processor. The @c value_type of
 * this property map will become the @c value_type of the distributed
 * property map. The distributed property map models the same property
 * map concepts as the @c LocalPropertyMap, with one exception: a
 * distributed property map cannot be an LvaluePropertyMap (because
 * remote values are not addressable), and is therefore limited to
 * ReadWritePropertyMap.
 */
template<typename GlobalMap, typename StorageMap,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen>,
         typename GhostCellStorage =
//                  detail::unordered_map_ghost_cell_storage<
                    detail::multi_index_ghost_cell_storage<
                      typename property_traits<GlobalMap>::key_type,
                      typename property_traits<StorageMap>::value_type>
        >
class distributed_property_map
{
 public:
  /// The key type of the property map.
  typedef typename property_traits<GlobalMap>::key_type key_type;

  /// The value type of the property map.
  typedef typename property_traits<StorageMap>::value_type value_type;
  typedef typename property_traits<StorageMap>::reference  reference;
  typedef typename amplusplus::transport::rank_type    rank_type;

 private:
  typedef distributed_property_map            self_type;
  typedef typename property_traits<StorageMap>::category local_category;
  typedef typename property_traits<StorageMap>::key_type local_key_type;
  typedef typename property_traits<GlobalMap>::value_type owner_local_pair;

  template <typename K, typename V, int N> friend class detail::associative_cache_ghost_cell_storage;
  template <typename K, typename V> friend class detail::multi_index_ghost_cell_storage;

 public:
  /** The property map category.  A distributed property map cannot be
   * an Lvalue Property Map, because values on remote processes cannot
   * be addresses.
   */
  typedef typename detail::make_nonlvalue_property_map<
    (is_base_and_derived<lvalue_property_map_tag, local_category>::value
     || is_same<lvalue_property_map_tag, local_category>::value)>
    ::template apply<local_category>::type category;

  /** Default-construct a distributed property map.  This function
   * creates an initialized property map that must be assigned to a
   * valid value before being used. It is only provided here because
   * property maps must be Default Constructible.
   */
  distributed_property_map() {}

  /** Construct a distributed property map.  Builds a distributed
   * property map communicating over the given process group and using
   * the given local property map for storage. Since no reduction
   * operation is provided, the default reduction operation @c
   * basic_reduce<value_type> is used.
   */
  distributed_property_map(amplusplus::transport& transport, const GlobalMap& global,
                           const StorageMap& pm, 
                           MessageGenerator message_gen = 
                             MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : data(new data_t(transport, global, pm, basic_reduce<value_type>(), 
                      basic_reduce<value_type>(), false, message_gen))
  {
    data->ghost_cells.reset(new GhostCellStorage());
    data->ghost_cells->unset_limit();
    data->setup_messages();
  }

  /** Construct a distributed property map.  Builds a distributed
   * property map communicating over the given process group and using
   * the given local property map for storage. The given @p reduce
   * parameter is used as the reduction operation.
   */
  template<typename Reduce>
  distributed_property_map(amplusplus::transport& transport, const GlobalMap& global,
                           const StorageMap& pm, const Reduce& reduce, 
                           MessageGenerator message_gen = 
                             MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)));

  ~distributed_property_map();

  /// Set the reduce operation of the distributed property map.
  template<typename Reduce>
  void set_reduce(const Reduce& reduce);

  // Set the consistency model for the distributed property map
  void set_consistency_model(int model);

  // Get the consistency model
  int get_consistency_model() const { return data->model; }

  // Set the maximum number of ghost cells that we are allowed to
  // maintain. If 0, all ghost cells will be retained.
  void set_max_ghost_cells(std::size_t max_ghost_cells);

  // Clear out all ghost cells
  void clear();

  // Reset the values in all ghost cells to the default value
  void reset();

  // Flush all values destined for remote processors
  void flush();

  reference operator[](const key_type& key) const
  {
    owner_local_pair p = get(data->global, key);
    
    if (p.first == data->transport.rank()) {
      return data->storage[p.second];
    } else {
      return cell(key);
    }
  }

  StorageMap&       base()       { return data->storage; }
  const StorageMap& base() const { return data->storage; }

  /** Sends a "put" request.
   * \internal
   *
   */
  void 
  request_put(rank_type, const key_type& k, const value_type& v) const
  { 
    data->put_msg.send(std::make_pair(k, v));
  }

  void 
  request_atomic_put(rank_type, const key_type& k, const value_type& v) const
  { 
    data->atomic_put_msg.send(std::make_pair(k, v));
  }

  /** Access the ghost cell for the given key.
   * \internal
   */
  value_type& cell(const key_type& k, bool request_if_missing = true) const;

  /** Perform synchronization
   * \internal
   */
  void synchronize();

  const GlobalMap& global() const { return data->global; }
  GlobalMap&       global()       { return data->global; }

  struct data_t
  {
    data_t(amplusplus::transport& trans, const GlobalMap& global, 
           const StorageMap& pm, const function1<value_type, key_type>& get_default_value,
           const function3<value_type, key_type, value_type, value_type>& reduce,
           bool has_default_resolver, MessageGenerator message_gen)
      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<std::pair<key_type, rank_type> >(),
          amplusplus::register_mpi_datatype<std::pair<key_type, value_type> >(),
          amplusplus::register_mpi_datatype<boost::tuple<rank_type, key_type, value_type> >(), 0)),
        transport(trans), global(global), storage(pm), 
        ghost_cells(), get_default_value(get_default_value), reduce(reduce), 
        has_default_resolver(has_default_resolver), model(cm_forward),
	put_msg(message_gen, trans, first_of_map_of_project1st<GlobalMap>(global), amplusplus::no_reduction),
	atomic_put_msg(message_gen, trans, first_of_map_of_project1st<GlobalMap>(global), amplusplus::no_reduction),
	get_msg(message_gen, trans, first_of_map_of_project1st<GlobalMap>(global), amplusplus::no_reduction),
	get_reply_msg(message_gen, trans, tuple_elem<0>(), amplusplus::no_reduction)
    { }

    const int dummy_first_member_for_init_order; // Unused

    /// The transport
    amplusplus::transport& transport;

    /// A mapping from the keys of this property map to the global
    /// descriptor.
    GlobalMap global;

    /// Local property map
    StorageMap storage;

    /// The ghost cells
    shared_ptr<GhostCellStorage> ghost_cells;

    /// Default value for remote ghost cells, as defined by the
    /// reduction operation.
    function1<value_type, key_type> get_default_value;

    /// Default value for remote ghost cells, as defined by the
    /// reduction operation.
    function3<value_type, key_type, value_type, value_type> reduce;

    /// True if this resolver is the "default" resolver, meaning that
    /// we should not be able to get() a default value; it needs to be
    /// request()ed first.
    bool has_default_resolver;

    // Current consistency model
    int model;

    // Function that resets all of the ghost cells to their default
    // values. It knows the type of the resolver, so we can eliminate
    // a large number of calls through function pointers.
    void (data_t::*reset)();

    // Clear out all ghost cells
    void clear();

    // Flush all values destined for remote processors
    void flush();

    // Send out requests to "refresh" the values of ghost cells that
    // we're holding.
    void refresh_ghost_cells();

  private:
    /// Incoming message handlers

    /** A request to store a value in a property map. The message
     * contains a std::pair<key, data>.
     */
    struct atomic_put_handler {
      
      explicit atomic_put_handler() : self(NULL) {}
      atomic_put_handler(data_t& self) : self(&self) {}
      
      void operator() (const std::pair<key_type, value_type>& req) const;
      
    protected:
      data_t* self;
    };

    /** A request to store a value in a property map. The message
     * contains a std::pair<key, data>.
     */
    struct put_handler {
      
      explicit put_handler() : self(NULL) {}
      put_handler(data_t& self) : self(&self) {}
      
      void operator() (const std::pair<key_type, value_type>& req) const;
      
    protected:
      data_t* self;
    };
    
    /** A request to retrieve a particular value in a property
     *  map. The message contains a key. The owner of that key will 
     *  reply with a value.
     */
    struct get_handler {
      
      explicit get_handler() : self(NULL) {}
      get_handler(data_t& self) : self(&self) {}
      
      void operator() (const std::pair<key_type, rank_type>& msg) const;
      
    protected:
      data_t* self;
    };
    
    struct get_reply_handler {
      
      explicit get_reply_handler() : self(NULL) {}
      get_reply_handler(data_t& self) : self(&self) {}
      
      void operator() (const boost::tuple<rank_type, key_type, value_type>& msg) const;
      
    protected:
      data_t* self;
    };
    
    /// Message types
    typedef typename MessageGenerator::template call_result<std::pair<key_type, value_type>, put_handler, 
							    first_of_map_of_project1st<GlobalMap>, 
							    amplusplus::no_reduction_t>::type
      put_message_type;
    // TODO: If reductions were defined at compile time we could perform local reductions here

    typedef typename MessageGenerator::template call_result<std::pair<key_type, value_type>, atomic_put_handler, 
							    first_of_map_of_project1st<GlobalMap>, 
							    amplusplus::no_reduction_t>::type
      atomic_put_message_type;

    typedef typename MessageGenerator::template call_result<std::pair<key_type, rank_type>, get_handler, 
							    first_of_map_of_project1st<GlobalMap>, 
							    amplusplus::no_reduction_t>::type
      get_message_type;

    typedef typename MessageGenerator::template call_result<boost::tuple<rank_type, key_type, value_type>, 
							    get_reply_handler, tuple_elem<0>, 
							    amplusplus::no_reduction_t>::type
      get_reply_message_type;

    /// Set handlers
    void setup_messages();
    
    /// Messages 
    put_message_type put_msg;
    atomic_put_message_type atomic_put_msg;  
    get_message_type get_msg;
    get_reply_message_type get_reply_msg;

  private:
    template<typename Resolver> void do_reset();

    friend class distributed_property_map;
  };
  friend struct data_t;

  shared_ptr<data_t> data;

 private:
  // Prunes the least recently used ghost cells until we have @c
  // max_ghost_cells or fewer ghost cells.
  void prune_ghost_cells() const;

  // Prunes a given ghost cell -- used by the ghost cell storage object.
  void prune_one_ghost_cell(const key_type&, const value_type&) const;
};

/* An implementation helper macro for the common case of naming
   distributed property maps with all of the normal template
   parameters. */
#define PBGL_DISTRIB_PMAP                                       \
  distributed_property_map<GlobalMap, StorageMap, MessageGenerator, GhostCellStorage>

/* The template declaration for a distributed property map. */
#define PBGL_DISTRIB_PMAP_TEMPLATE \
  template<typename GlobalMap, typename StorageMap, typename MessageGenerator, \
           typename GhostCellStorage>

/* Request that the value for the given remote key be retrieved in
   the next synchronization round. */
PBGL_DISTRIB_PMAP_TEMPLATE
inline void
request(const PBGL_DISTRIB_PMAP& pm,
        typename PBGL_DISTRIB_PMAP::key_type const& key)
{
  if (get(pm.data->global, key).first != pm.data->transport.rank())
    pm.cell(key, true);
}

/** Get the value associated with a particular key.  Retrieves the
 * value associated with the given key. If the key denotes a
 * locally-owned object, it returns the value from the local property
 * map; if the key denotes a remotely-owned object, retrieves the
 * value of the ghost cell for that key, which may be the default
 * value provided by the reduce operation.
 *
 * Complexity: For a local key, O(1) get operations on the underlying
 * property map. For a non-local key, O(1) accesses to the ghost cells.
 */
PBGL_DISTRIB_PMAP_TEMPLATE
inline
typename PBGL_DISTRIB_PMAP::value_type
get(const PBGL_DISTRIB_PMAP& pm,
    typename PBGL_DISTRIB_PMAP::key_type const& key)
{
  using boost::get;

  typename property_traits<GlobalMap>::value_type p = 
    get(pm.data->global, key);

  if (p.first == pm.data->transport.rank()) {
    return get(pm.data->storage, p.second);
  } else {
    return pm.cell(key);
  }
}

/** Put a value associated with the given key into the property map.
 * When the key denotes a locally-owned object, this operation updates
 * the underlying local property map. Otherwise, the local ghost cell
 * is updated and a "put" message is sent to the processor owning this
 * key.
 *
 * Complexity: For a local key, O(1) put operations on the underlying
 * property map. For a nonlocal key, O(1) accesses to the ghost cells
 * and will send O(1) messages of size O(sizeof(key) + sizeof(value)).
 */
PBGL_DISTRIB_PMAP_TEMPLATE
void
put(const PBGL_DISTRIB_PMAP& pm,
    typename PBGL_DISTRIB_PMAP::key_type const & key,
    typename PBGL_DISTRIB_PMAP::value_type const & value)
{
  using boost::put;

  typename property_traits<GlobalMap>::value_type p = 
    get(pm.data->global, key);

  if (p.first == pm.data->transport.rank()) {
    put(pm.data->storage, p.second, value);
    // NGE: Apply reduce?
//     detail::maybe_put(pm.data->storage, p.second, 
// 		      pm.data->reduce(key,
// 				      get(pm.data->storage, p.second),
// 				      value));
  } else {
    if (pm.data->model & cm_forward) 
      pm.request_put(p.first, key, value);
#ifndef THREADSAFE_GHOST_CELLS
    // NGE: If there is a reduction operation we might need to make
    //      sure cell gets the correct value.  Unfortunately we have 
    //      know way of knowing whether this caused a new cell to be 
    //      inserted or not.  If it's a new cell then the value returned
    //      is the value given by the default resolver, otherwise it's
    //      the existing ghost cell value.
    typename PBGL_DISTRIB_PMAP::value_type& new_cell = pm.cell(key, false);
    new_cell = pm.data->reduce(key, new_cell, value);
#endif
  }
}

/** Put a value associated with the given key into the property map.
 * When the key denotes a locally-owned object, this operation updates
 * the underlying local property atomically map. Otherwise, an atomic
 * "put" message is sent to the processor owning this key.
 *
 * Complexity: For a local key, O(infinity) exchange operations on the 
 * underlying property map. For a nonlocal key, O(1) messages of size 
 * O(sizeof(key) + sizeof(value)).
 */
PBGL_DISTRIB_PMAP_TEMPLATE
void
atomic_put(const PBGL_DISTRIB_PMAP& pm,
           typename PBGL_DISTRIB_PMAP::key_type const & key,
           typename PBGL_DISTRIB_PMAP::value_type const & value)
{
  using boost::put;

  typename property_traits<GlobalMap>::value_type p = 
    get(pm.data->global, key);

  if (p.first == pm.data->transport.rank()) {

    typename PBGL_DISTRIB_PMAP::value_type oldval, newval;

    do {
      oldval = get(pm.data->storage, p.second);
      newval = pm.data->reduce(key, oldval, value);
    } while(!maybe_exchange(pm.data->storage, p.second, oldval, newval));
  } else {
    if (pm.data->model & cm_forward) 
      pm.request_atomic_put(p.first, key, value);
#ifndef THREADSAFE_GHOST_CELLS
    pm.cell(key, false) = value;
#endif
  }
}

/** Put a value associated with a given key into the local view of the
 * property map. This operation is equivalent to @c put, but with one
 * exception: no message will be sent to the owning processor in the
 * case of a remote update. The effect is that any value written via
 * @c local_put for a remote key may be overwritten in the next
 * synchronization round.
 */
PBGL_DISTRIB_PMAP_TEMPLATE
void
local_put(const PBGL_DISTRIB_PMAP& pm,
          typename PBGL_DISTRIB_PMAP::key_type const & key,
          typename PBGL_DISTRIB_PMAP::value_type const & value)
{
  using boost::put;

  typename property_traits<GlobalMap>::value_type p = 
    get(pm.data->global, key);

  if (p.first == pm.data->transport.rank())
    put(pm.data->storage, p.second, value);
  else pm.cell(key, false) = value;
}

/** Cache the value associated with the given remote key. If the key
 *  is local, ignore the operation. */
PBGL_DISTRIB_PMAP_TEMPLATE
inline void
cache(const PBGL_DISTRIB_PMAP& pm,
      typename PBGL_DISTRIB_PMAP::key_type const & key,
      typename PBGL_DISTRIB_PMAP::value_type const & value)
{
  typename PBGL_DISTRIB_PMAP::rank_type id = get(pm.data->global, key).first;

  if (id != pm.data->transport.rank()) pm.cell(key, false) = value;
}

/** Exchange (i.e. Compare-and-swap) the old value of the given local
 *  key with the new value. */
PBGL_DISTRIB_PMAP_TEMPLATE
inline bool
exchange(const PBGL_DISTRIB_PMAP& pm,
         typename PBGL_DISTRIB_PMAP::key_type const & key,
         typename PBGL_DISTRIB_PMAP::value_type const & oldval,
         typename PBGL_DISTRIB_PMAP::value_type const & newval)
{
  // Exchange is only valid on local data until we get better one-sided ops
  assert(get(pm.data->global, key).first == pm.data->transport.rank()); 

  return exchange(pm.data->storage, get(pm.data->global, key).second, oldval, newval);
}

/** Exchange (i.e. Compare-and-swap) the old value of the given local
 *  key with the new value if CAS is supported. */
PBGL_DISTRIB_PMAP_TEMPLATE
inline bool
maybe_exchange(const PBGL_DISTRIB_PMAP& pm,
               typename PBGL_DISTRIB_PMAP::key_type const & key,
               typename PBGL_DISTRIB_PMAP::value_type const & oldval,
               typename PBGL_DISTRIB_PMAP::value_type const & newval)
{
  // Exchange is only valid on local data until we get better one-sided ops
  assert(get(pm.data->global, key).first == pm.data->transport.rank()); 

  return maybe_exchange(pm.data->storage, get(pm.data->global, key).second, oldval, newval);
}

/// Synchronize the property map.
/// Must never be called *during* an epoch
PBGL_DISTRIB_PMAP_TEMPLATE
void
synchronize(PBGL_DISTRIB_PMAP& pm)
{
  pm.synchronize();
}

/// Create a distributed property map.
PBGL_DISTRIB_PMAP_TEMPLATE
inline PBGL_DISTRIB_PMAP
make_distributed_property_map(amplusplus::transport& trans, GlobalMap global, 
                              StorageMap storage)
{
  typedef PBGL_DISTRIB_PMAP result_type;
  return result_type(trans, global, storage);
}

/**
 * \overload
 */
template<typename GlobalMap, typename StorageMap, typename MessageGenerator,
         typename GhostCellStorage, typename Reduce>
inline PBGL_DISTRIB_PMAP
make_distributed_property_map(amplusplus::transport& trans, GlobalMap global, 
                              StorageMap storage, Reduce reduce)
{
  typedef PBGL_DISTRIB_PMAP result_type;
  return result_type(trans, global, storage, reduce);
}

} } // end namespace boost::parallel

// Boost's functional/hash
namespace boost {
  template<typename D, typename V>
  struct hash<boost::detail::edge_desc_impl<D, V> >
  {
    std::size_t operator()(const boost::detail::edge_desc_impl<D, V> & x) const
    { return hash_value(x.get_property()); }
  };
}

#include <boost/property_map/parallel/impl/distributed_property_map.cpp>

#undef PBGL_DISTRIB_PMAP
#undef PBGL_DISTRIB_PMAP_TEMPLATE

/**
 *  All of this used to live at the end of property_map.hpp and
 *  vector_property_map.hpp but since AM++ uses property maps
 *  internally it had to be moved to a place that could be included
 *  *after* AM++ headers.
 */

#include <boost/property_map/parallel/local_property_map.hpp>

namespace boost {

  /** Exchange method for property maps
   * 
   * Add exchange() method to sequential property maps
   *
   */
  template <class PropertyMap, class Reference, class K, class V>
  inline bool
  exchange(const put_get_helper<Reference, PropertyMap>& pa, K k, 
           const V& oldval, const V& newval)
  {
    using boost::parallel::bool_compare_and_swap;

    return bool_compare_and_swap(&static_cast<const PropertyMap&>(pa)[k], 
                                 oldval, newval);
  }

  // 
  // Use CAS if available, otherwise put, locking will be handled by the caller
  //
  template <class PropertyMap, class Reference, class K, class V>
  inline typename boost::enable_if<boost::parallel::atomics_supported<V>, bool>::type
  maybe_exchange(const put_get_helper<Reference, PropertyMap>& pa, K k, 
                 const V& oldval, const V& newval)
  {
    return exchange(pa, k, oldval, newval);
  } 
    
  template <class PropertyMap, class Reference, class K, class V>
  inline typename boost::disable_if<boost::parallel::atomics_supported<V>, bool>::type
  maybe_exchange(const put_get_helper<Reference, PropertyMap>& pa, K k, 
                 const V& oldval, const V& newval)
  {
    if (static_cast<const PropertyMap&>(pa)[k] == oldval) {
      static_cast<const PropertyMap&>(pa)[k] = newval;
      return true;
    }
    
    return false;
  } 

/** Distributed iterator property map.
 *
 * This specialization of @ref iterator_property_map builds a
 * distributed iterator property map given the local index maps
 * generated by distributed graph types that automatically have index
 * properties. 
 *
 * This specialization is useful when creating external distributed
 * property maps via the same syntax used to create external
 * sequential property maps.
 */
template<typename RandomAccessIterator, typename GlobalMap, 
         typename StorageMap, typename ValueType, typename Reference>
class iterator_property_map
        <RandomAccessIterator, 
         local_property_map<GlobalMap, StorageMap>,
         ValueType, Reference>
  : public parallel::distributed_property_map
             <GlobalMap, 
              iterator_property_map<RandomAccessIterator, StorageMap,
                                    ValueType, Reference> >
{
  typedef iterator_property_map<RandomAccessIterator, StorageMap, 
                                ValueType, Reference> local_iterator_map;

  typedef parallel::distributed_property_map<GlobalMap, local_iterator_map> inherited;

  typedef local_property_map<GlobalMap, StorageMap> index_map_type;
  typedef iterator_property_map self_type;

public:
  iterator_property_map() { }

  iterator_property_map(RandomAccessIterator cc, const index_map_type& id)
    : inherited(id.transport(), id.global(), 
                local_iterator_map(cc, id.base())) { }
};

/** Distributed iterator property map.
 *
 * This specialization of @ref iterator_property_map builds a
 * distributed iterator property map given a distributed index
 * map. Only the local portion of the distributed index property map
 * is utilized.
 *
 * This specialization is useful when creating external distributed
 * property maps via the same syntax used to create external
 * sequential property maps.
 */
template<typename RandomAccessIterator, typename GlobalMap, 
         typename StorageMap, typename ValueType, typename Reference>
class iterator_property_map<
        RandomAccessIterator, 
        parallel::distributed_property_map<GlobalMap, StorageMap>,
        ValueType, Reference
      >
  : public parallel::distributed_property_map
             <GlobalMap,
              iterator_property_map<RandomAccessIterator, StorageMap,
                                    ValueType, Reference> >
{
  typedef iterator_property_map<RandomAccessIterator, StorageMap,
                                ValueType, Reference> local_iterator_map;

  typedef parallel::distributed_property_map<GlobalMap, local_iterator_map> inherited;

  typedef parallel::distributed_property_map<GlobalMap, StorageMap> index_map_type;

public:
  iterator_property_map() { }

  iterator_property_map(RandomAccessIterator cc, const index_map_type& id)
    : inherited(id.transport(), id.global(),
                local_iterator_map(cc, id.base())) { }
};

namespace parallel {
// Generate an iterator property map with a specific kind of ghost
// cells
template<typename RandomAccessIterator, typename GlobalMap, typename StorageMap>
distributed_property_map<GlobalMap,
                         iterator_property_map<RandomAccessIterator, 
                                               StorageMap> >
make_iterator_property_map(RandomAccessIterator cc,
                           local_property_map<GlobalMap, StorageMap> index_map)
{
  typedef distributed_property_map<GlobalMap,
            iterator_property_map<RandomAccessIterator, StorageMap> >
    result_type;
  return result_type(index_map.transport(), index_map.global(),
                     make_iterator_property_map(cc, index_map.base()));
}

} // end namespace parallel

/** Distributed safe iterator property map.
 *
 * This specialization of @ref safe_iterator_property_map builds a
 * distributed iterator property map given the local index maps
 * generated by distributed graph types that automatically have index
 * properties. 
 *
 * This specialization is useful when creating external distributed
 * property maps via the same syntax used to create external
 * sequential property maps.
 */
template<typename RandomAccessIterator, typename GlobalMap, 
         typename StorageMap, typename ValueType, typename Reference>
class safe_iterator_property_map
        <RandomAccessIterator, 
         local_property_map<GlobalMap, StorageMap>,
         ValueType, Reference>
  : public parallel::distributed_property_map
             <GlobalMap,
              safe_iterator_property_map<RandomAccessIterator, StorageMap,
                                         ValueType, Reference> >
{
  typedef safe_iterator_property_map<RandomAccessIterator, StorageMap, 
                                     ValueType, Reference> local_iterator_map;

  typedef parallel::distributed_property_map<GlobalMap, local_iterator_map> inherited;

  typedef local_property_map<GlobalMap, StorageMap> index_map_type;

public:
  safe_iterator_property_map() { }

  safe_iterator_property_map(RandomAccessIterator cc, std::size_t n, 
                             const index_map_type& id, int coalescing_size = 1 << 12)
    : inherited(id.transport(), id.global(),
                local_iterator_map(cc, n, id.base()), coalescing_size) { }
};

/** Distributed safe iterator property map.
 *
 * This specialization of @ref safe_iterator_property_map builds a
 * distributed iterator property map given a distributed index
 * map. Only the local portion of the distributed index property map
 * is utilized.
 *
 * This specialization is useful when creating external distributed
 * property maps via the same syntax used to create external
 * sequential property maps.
 */
template<typename RandomAccessIterator, typename GlobalMap, 
         typename StorageMap, typename ValueType, typename Reference>
class safe_iterator_property_map<
        RandomAccessIterator, 
        parallel::distributed_property_map<GlobalMap, StorageMap>,
        ValueType, Reference>
  : public parallel::distributed_property_map
             <GlobalMap,
              safe_iterator_property_map<RandomAccessIterator, StorageMap,
                                         ValueType, Reference> >
{
  typedef safe_iterator_property_map<RandomAccessIterator, StorageMap,
                                     ValueType, Reference> local_iterator_map;

  typedef parallel::distributed_property_map<GlobalMap, local_iterator_map> inherited;

  typedef parallel::distributed_property_map<GlobalMap, StorageMap> index_map_type;

public:
  safe_iterator_property_map() { }

  safe_iterator_property_map(RandomAccessIterator cc, std::size_t n, 
                             const index_map_type& id, int coalescing_size = 1 << 12)
    : inherited(id.transport(), id.global(), 
                local_iterator_map(cc, n, id.base()), coalescing_size) { }
};                                            

/** Distributed vector property map.
 *
 * This specialization of @ref vector_property_map builds a
 * distributed vector property map given the local index maps
 * generated by distributed graph types that automatically have index
 * properties. 
 *
 * This specialization is useful when creating external distributed
 * property maps via the same syntax used to create external
 * sequential property maps.
 */
template<typename T, typename GlobalMap, typename StorageMap>
class vector_property_map<T, local_property_map<GlobalMap, StorageMap> >
  : public parallel::distributed_property_map<
             GlobalMap, vector_property_map<T, StorageMap> >
{
  typedef vector_property_map<T, StorageMap> local_iterator_map;

  typedef parallel::distributed_property_map<GlobalMap, local_iterator_map> inherited;

  typedef local_property_map<GlobalMap, StorageMap> index_map_type;

public:
  vector_property_map(const index_map_type& index = index_map_type())
    : inherited(index.transport(), index.global(),
                local_iterator_map(index.base())) { }

  vector_property_map(unsigned inital_size, 
                      const index_map_type& index = index_map_type())
    : inherited(index.transport(),  index.global(),
                local_iterator_map(inital_size, index.base())) { }
};

/** Distributed vector property map.
 *
 * This specialization of @ref vector_property_map builds a
 * distributed vector property map given the local index maps
 * generated by distributed graph types that automatically have index
 * properties. 
 *
 * This specialization is useful when creating external distributed
 * property maps via the same syntax used to create external
 * sequential property maps.
 */
template<typename T, typename GlobalMap, typename StorageMap>
class vector_property_map<
        T, 
        parallel::distributed_property_map<
          GlobalMap,
          StorageMap
        >
      > 
  : public parallel::distributed_property_map<
             GlobalMap, vector_property_map<T, StorageMap> >
{
  typedef vector_property_map<T, StorageMap> local_iterator_map;

  typedef parallel::distributed_property_map<GlobalMap, local_iterator_map> 
    inherited;

  typedef parallel::distributed_property_map<GlobalMap, StorageMap>
    index_map_type;

public:
  vector_property_map(const index_map_type& index = index_map_type())
    : inherited(index.transport(), index.global(),
                local_iterator_map(index.base())) { }

  vector_property_map(unsigned inital_size, 
                      const index_map_type& index = index_map_type())
    : inherited(index.transport(), index.global(),
                local_iterator_map(inital_size, index.base())) { }
};

// Import distributed property map put/get for vector and iterator
// property maps which inherit from distributed property map
using boost::parallel::put;
using boost::parallel::get;

} // namespace boost

// These need to go in global namespace because Koenig
// lookup does not apply to T*.

// V must be convertible to T
template <class T, class V>
inline bool exchange(T* pa, std::ptrdiff_t k, const V& oldval, const V& newval) 
{ 
  using boost::parallel::bool_compare_and_swap;
  
  return bool_compare_and_swap(&pa[k], oldval, newval);
}

// 
// Use CAS if available, otherwise put, locking will be handled by the caller
//
template <class T, class V>
inline typename boost::enable_if<boost::parallel::atomics_supported<T>, bool>::type
maybe_exchange(T* pa, std::ptrdiff_t k, const V& oldval, const V& newval)
{
  return exchange(pa, k, oldval, newval);
} 

template <class T, class V>
inline typename boost::disable_if<boost::parallel::atomics_supported<T>, bool>::type
maybe_exchange(T* pa, std::ptrdiff_t k, const V& oldval, const V& newval)
{
  if (get(pa, k) == oldval) {
    put(pa, k, newval);
    return true;
  }
  
  return false;
} 

#endif // BOOST_PARALLEL_DISTRIBUTED_PROPERTY_MAP_HPP
