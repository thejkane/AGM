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
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/test/minimal.hpp>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/graph/parallel/basic_reduce.hpp>
#include <boost/functional/hash.hpp>

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

using namespace boost;

// enum color_t { red, blue };

typedef unsigned char color_t;
color_t red = 0x01;
color_t blue = 0x02;

struct remote_key 
{
  remote_key(std::size_t p = std::numeric_limits<std::size_t>::max(), std::size_t l = 0) 
    : processor(p), local_key(l) {}

  std::size_t processor;
  std::size_t local_key;
};

namespace amplusplus {
  template<>
  struct make_mpi_datatype<remote_key> : make_mpi_datatype_base {
    make_mpi_datatype<std::size_t> dt1;
    make_mpi_datatype<std::size_t> dt2;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1(), dt2() {
      int blocklengths[2] = {1, 1};
      MPI_Aint displacements[2];
      remote_key test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.processor, &displacements[0]);
      MPI_Get_address(&test_object.local_key, &displacements[1]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      MPI_Datatype types[2] = {dt1.get(), dt2.get()};
      MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
}

namespace boost { 

template<>
struct hash<remote_key>
{
  std::size_t operator()(const remote_key& key) const
  {
    std::size_t hash = hash_value(key.processor);
    hash_combine(hash, key.local_key);
    return hash;
  }
};

}

inline bool operator==(const remote_key& x, const remote_key& y)
{ return x.processor == y.processor && x.local_key == y.local_key; }

struct remote_key_to_global
{
  typedef readable_property_map_tag category;
  typedef remote_key key_type;
  typedef std::pair<std::size_t, std::size_t> value_type;
  typedef value_type reference;
};

inline std::pair<std::size_t, std::size_t> 
get(remote_key_to_global, const remote_key& key)
{
  return std::make_pair(key.processor, key.local_key);
}

template<typename T>
struct my_reduce : boost::parallel::basic_reduce<T> {
  BOOST_STATIC_CONSTANT(bool, non_default_resolver = true);
};

void colored_test(amplusplus::transport& trans)
{
  const unsigned int n = 500;
  
  color_t my_start_color = trans.rank() % 2 == 0? ::red : ::blue;
  unsigned int next_processor = (trans.rank() + 1) % trans.size();
  color_t next_start_color = next_processor % 2 == 0? ::red : ::blue;

  // Initial color map: even-numbered processes are all red,
  // odd-numbered processes are all blue.
  std::vector<color_t> color_vec(n, my_start_color);

  typedef iterator_property_map<std::vector<color_t>::iterator, 
                                identity_property_map> LocalPropertyMap;
  LocalPropertyMap local_colors(color_vec.begin(), identity_property_map());

  // Create the distributed property map
  typedef boost::parallel::distributed_property_map<remote_key_to_global,
                                                    LocalPropertyMap> ColorMap;
  ColorMap colors(trans, remote_key_to_global(), local_colors);
  colors.set_reduce(my_reduce<color_t>());

  if (trans.rank() == 0) std::cerr << "Checking local colors...";
  // check local processor colors
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(trans.rank(), i);
    BOOST_CHECK(get(colors, k) == my_start_color);
  }

  colors.set_consistency_model(boost::parallel::cm_bidirectional);
  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking next processor's default colors...";
  // check next processor's colors
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(colors, k) == color_t());
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking next processor's colors...";
  // Check next pocessor's colors when requests complete
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(colors, k) == next_start_color); // colors were requested above
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChanging next processor's colors...";
  // change the next processor's colors
  color_t next_finish_color = next_processor % 2 == 0? ::blue : ::red;
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    put(colors, k, next_finish_color);
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();


  if (trans.rank() == 0) std::cerr << "OK.\nChecking changed colors...";
  // check our own colors
  color_t my_finish_color = trans.rank() % 2 == 0? ::blue : ::red;
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(trans.rank(), i);
    BOOST_CHECK(get(colors, k) == my_finish_color);
  }

  synchronize(colors); // refresh ghost cells
  trans.begin_epoch();

  // check our neighbor's colors
  if (trans.rank() == 0) std::cerr << "OK.\nChecking changed colors on neighbor...";
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(colors, k) == next_finish_color);
  }

  trans.end_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\n";
}

#if 0 // no datatype for bools and std:string is variable sized which requires serialization
void bool_test(amplusplus::transport& trans)
{
  const unsigned int n = 500;
  
  bool my_start_value = trans.rank() % 2;
  unsigned int next_processor = (trans.rank() + 1) % trans.size();
  bool next_start_value = !my_start_value;

  // Initial color map: even-numbered processes are false, 
  // odd-numbered processes are true
  std::vector<bool> bool_vec(n, my_start_value);

  typedef iterator_property_map<std::vector<bool>::iterator, 
                                identity_property_map> LocalPropertyMap;
  LocalPropertyMap local_values(bool_vec.begin(), identity_property_map());

  // Create the distributed property map
  typedef boost::parallel::distributed_property_map<remote_key_to_global,
                                                    LocalPropertyMap> ValueMap;
  ValueMap values(trans, remote_key_to_global(), local_values);
  values.set_reduce(my_reduce<bool>());

  if (trans.rank() == 0) std::cerr << "Checking local values...";
  // check local processor values
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(trans.rank(), i);
    BOOST_CHECK(get(values, k) == my_start_value);
  }

  values.set_consistency_model(boost::parallel::cm_bidirectional);
  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking next processor's default values...";
  // check next processor's values
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(values, k) == false);
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking next processor's values...";
  // check next processor's values
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(values, k) == next_start_value); // values were requested in previous loop
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChanging next processor's values...";
  // change the next processor's values
  bool next_finish_value = next_processor % 2 == 0;
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    put(values, k, next_finish_value);
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking changed values...";
  // check our own values
  bool my_finish_value = trans.rank() % 2 == 0;
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(trans.rank(), i);
    BOOST_CHECK(get(values, k) == my_finish_value);
  }

  trans.begin_epoch();

  synchronize(values); // refresh ghost cells
  trans.begin_epoch();

  // check our neighbor's values
  if (trans.rank() == 0) std::cerr << "OK.\nChecking changed values on neighbor...";
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(values, k) == next_finish_value);
  }

  trans.end_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\n";
}

template <typename Transport>
void string_test(Transport& trans)
{
  const unsigned int n = 500;
  
  std::string my_start_string = lexical_cast<std::string>(trans.rank());
  unsigned int next_processor = (trans.rank() + 1) % trans.size();
  std::string next_start_string = lexical_cast<std::string>(next_processor);

  // Initial color map: even-numbered processes are false, 
  // odd-numbered processes are true
  std::vector<std::string> string_vec(n, my_start_string);

  typedef iterator_property_map<std::vector<std::string>::iterator, 
                                identity_property_map> LocalPropertyMap;
  LocalPropertyMap local_strings(string_vec.begin(), identity_property_map());

  // Create the distributed property map
  typedef boost::parallel::distributed_property_map<Transport,
                                                    remote_key_to_global,
                                                    LocalPropertyMap> StringMap;
  StringMap strings(trans, remote_key_to_global(), local_strings);
  strings.set_reduce(my_reduce<std::string>());

  if (trans.rank() == 0) std::cerr << "Checking local strings...";
  // check local processor strings
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(trans.rank(), i);
    BOOST_CHECK(get(strings, k) == my_start_string);
  }

  strings.set_consistency_model(boost::parallel::cm_bidirectional);
  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking next processor's default strings...";
  // check next processor's strings
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(strings, k) == std::string());
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking next processor's strings...";
  // check next processor's strings
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(strings, k) == next_start_string);
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  trans.begin_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChanging next processor's strings...";
  // change the next processor's strings
  std::string next_finish_string = next_start_string + next_start_string;
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    put(strings, k, next_finish_string);
  }

  if (trans.rank() == 0) std::cerr << "OK.\nSynchronizing...";
  trans.end_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\nChecking changed strings...";
  // check our own strings
  std::string my_finish_string = my_start_string + my_start_string;
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(trans.rank(), i);
    BOOST_CHECK(get(strings, k) == my_finish_string);
  }

  synchronize(strings); // refresh ghost cells
  trans.begin_epoch();

  // check our neighbor's strings
  if (trans.rank() == 0) std::cerr << "OK.\nChecking changed strings on neighbor...";
  for (unsigned int i = 0; i < n; ++i) {
    remote_key k(next_processor, i);
    BOOST_CHECK(get(strings, k) == next_finish_string);
  }

  trans.end_epoch();

  if (trans.rank() == 0) std::cerr << "OK.\n";
}
#endif

int test_main(int argc, char** argv)
{
  amplusplus::environment env = amplusplus::mpi_environment(argc, argv);
  amplusplus::transport trans = env.create_transport();

  trans.set_nthreads(1);

  amplusplus::register_mpi_datatype<remote_key>();

  colored_test(trans);
  // bool_test(trans);
  // string_test(trans);
  return 0;
}
