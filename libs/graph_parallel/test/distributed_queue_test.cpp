// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nick Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>

#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/distributed/queue.hpp>
#include <boost/test/minimal.hpp>
#include <boost/pending/queue.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <utility>
#include <iostream>

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

struct global_value 
{
  global_value(int p = -1, std::size_t l = 0) : processor(p), value(l) {}

  int processor;
  std::size_t value;
};

inline bool operator==(const global_value& x, const global_value& y)
{ return x.processor == y.processor && x.value == y.value; }

struct global_value_owner_map
{
  typedef int value_type;
  typedef value_type reference;
  typedef global_value key_type;
  typedef boost::readable_property_map_tag category;
};

global_value_owner_map::value_type get(global_value_owner_map, global_value k)
{
  return k.processor;
}

namespace amplusplus {
  // global_value 
  template <>
  struct make_mpi_datatype<global_value> : make_mpi_datatype_base {
    make_mpi_datatype<int> dt1;
    make_mpi_datatype<std::size_t> dt2;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1(), dt2() {
      int blocklengths[2] = {1, 1};
      MPI_Aint displacements[2];
      global_value test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.processor, &displacements[0]);
      MPI_Get_address(&test_object.value, &displacements[1]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      MPI_Datatype types[2] = {dt1.get(), dt2.get()};
      MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}

void test_distributed_queue(amplusplus::transport& trans)
{
  amplusplus::register_mpi_datatype<global_value>();

  typedef append_buffer<global_value> local_queue_type;

  typedef boost::graph::distributed::distributed_queue<global_value_owner_map,
                                                       local_queue_type> dist_queue_type;

  amplusplus::transport::rank_type id = trans.rank(), n = trans.size();

  //
  // Split-phase queue test
  //

  dist_queue_type Q(trans, global_value_owner_map(), 
		    boost::make_shared<append_buffer<global_value> >(),
		    boost::make_shared<append_buffer<global_value> >());

  global_value v(0, 0);

  if (id == 0) {
    std::cerr << "Should print level of each processor in a binary tree:\n";
  }

  {
    amplusplus::scoped_epoch epoch(trans);
    if (id == n-1) Q.push(v);
  }

  unsigned long all_queues_empty;

  while (true) {
    Q.swap();
    dist_queue_type::size_type starting_size = Q.size();

    {
      amplusplus::scoped_epoch_value epoch(trans, (starting_size == 0), all_queues_empty);

      if (starting_size) {
// 	v = Q.top(); Q.pop();
	v = Q[0]; Q.clear();
	
	std::cerr << "#" << id << ": level = " << v.value << std::endl;
	
	std::size_t level_begin = 1;
	for (std::size_t i = 0; i < v.value; ++i) level_begin *= 2;
	std::size_t level_end = level_begin * 2;
	BOOST_CHECK(level_begin <= (id + 1));
	BOOST_CHECK((id + 1) <= level_end);
	
	++v.value;
	v.processor = v.processor * 2 + 1;
	
	if (v.processor < (int)n) Q.push(v);
	++v.processor;
	if (v.processor < (int)n) Q.push(v);
      }
    }
    if (all_queues_empty == trans.size()) break;
  }

  // TODO (NGE): Non-split phase queue test, need appropriate local queue type
}

int test_main(int argc, char** argv)
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(1);
  test_distributed_queue(trans);
  return 0;
}
