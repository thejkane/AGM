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

#ifndef AMPLUSPLUS_MAKE_MPI_DATATYPE_HPP
#define AMPLUSPLUS_MAKE_MPI_DATATYPE_HPP

#include <mpi.h>
#include <utility>
#include <boost/config.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/noncopyable.hpp>
#include <am++/detail/type_info_map.hpp>
#ifndef BOOST_NO_0X_HDR_TUPLE
#include <tuple>
#endif

namespace amplusplus {

class scoped_mpi_datatype: boost::noncopyable {
  MPI_Datatype dt;
  public:
  explicit scoped_mpi_datatype(): dt(MPI_DATATYPE_NULL) {}
  // explicit scoped_mpi_datatype(MPI_Datatype dt = MPI_DATATYPE_NULL): dt(dt) {}
  operator MPI_Datatype() const {return dt;}
  MPI_Datatype* get_ptr() {return &dt;}
  MPI_Datatype get() const {return dt;}
  ~scoped_mpi_datatype() {if (dt != MPI_DATATYPE_NULL) MPI_Type_free(&dt);}
};

struct make_mpi_datatype_base {
  virtual ~make_mpi_datatype_base() {}
  virtual MPI_Datatype get() const = 0;
};

template <typename T, typename Enable = void>
struct make_mpi_datatype {};

template <> struct make_mpi_datatype<char> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_CHAR;}};
template <> struct make_mpi_datatype<signed char> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_CHAR;}};
template <> struct make_mpi_datatype<unsigned char> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_UNSIGNED_CHAR;}};
template <> struct make_mpi_datatype<short> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_SHORT;}};
template <> struct make_mpi_datatype<unsigned short> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_UNSIGNED_SHORT;}};
template <> struct make_mpi_datatype<int> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_INT;}};
template <> struct make_mpi_datatype<unsigned int> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_UNSIGNED;}};
template <> struct make_mpi_datatype<long> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_LONG;}};
template <> struct make_mpi_datatype<unsigned long> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_UNSIGNED_LONG;}};
template <> struct make_mpi_datatype<long long> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_LONG_LONG;}};
template <> struct make_mpi_datatype<unsigned long long> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_UNSIGNED_LONG_LONG;}};
template <> struct make_mpi_datatype<float> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_FLOAT;}};
template <> struct make_mpi_datatype<double> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_DOUBLE;}};
template <> struct make_mpi_datatype<long double> : make_mpi_datatype_base {MPI_Datatype get() const {return MPI_LONG_DOUBLE;}};

template <typename A, typename B>
struct make_mpi_datatype<std::pair<A, B> > : make_mpi_datatype_base {
  make_mpi_datatype<A> dt1;
  make_mpi_datatype<B> dt2;
  scoped_mpi_datatype dt;
  make_mpi_datatype(): dt1(), dt2() {
    int blocklengths[2] = {1, 1};
    MPI_Aint displacements[2];
    char dummy;
    std::pair<A, B> *test_object = (std::pair<A, B>*)(&dummy);
    MPI_Aint test_object_ptr;
    MPI_Get_address(test_object, &test_object_ptr);
    MPI_Get_address(&test_object->first, &displacements[0]);
    MPI_Get_address(&test_object->second, &displacements[1]);
    displacements[0] -= test_object_ptr;
    displacements[1] -= test_object_ptr;
    MPI_Datatype types[2] = {dt1.get(), dt2.get()};
    MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
    MPI_Type_commit(dt.get_ptr());
  }
  MPI_Datatype get() const {return dt;}
};

namespace detail {
  namespace tuple_check {
    BOOST_MPL_HAS_XXX_TRAIT_DEF(inherited);

    template <typename T>
    struct is_actual_cons: public boost::mpl::false_ {};

    template <typename Hd, typename Tl>
    struct is_actual_cons<boost::tuples::cons<Hd, Tl> >
      : public boost::mpl::true_ {};

    template <typename T>
    struct inherited_is_cons {
      typedef boost::mpl::bool_<is_actual_cons<typename T::inherited>::value> type;
    };

    template <typename T>
    struct is_cons
      : public boost::mpl::eval_if<
                 has_inherited<T>,
                 inherited_is_cons<T>,
                 boost::mpl::false_>
      {};

    template <typename T>
    struct is_nil
      : public boost::is_base_and_derived<boost::tuples::null_type, T> {};
  }

  template <typename T>
  struct is_boost_tuple
    : public boost::mpl::or_<tuple_check::is_cons<T>,
                             tuple_check::is_nil<T> >
    {};
}

template <typename T>
struct make_mpi_datatype<T,
                         typename boost::enable_if<
                                    detail::is_boost_tuple<T>
                                  >::type>
       : make_mpi_datatype_base
{
  make_mpi_datatype<typename T::inherited> dt_inh;
  make_mpi_datatype(): dt_inh() {}
  MPI_Datatype get() const {return dt_inh.get();}
};

template <typename Hd, typename Tl>
struct make_mpi_datatype<boost::tuples::cons<Hd, Tl> > : make_mpi_datatype_base {
  make_mpi_datatype<Hd> dt_hd;
  make_mpi_datatype<Tl> dt_tl;
  scoped_mpi_datatype dt;

  make_mpi_datatype(): dt_hd(), dt_tl() {
    int blocklengths[2] = {1, 1};
    MPI_Aint displacements[2];
    char dummy;
    boost::tuples::cons<Hd, Tl> *test_object = (boost::tuples::cons<Hd, Tl>*)(&dummy);
    MPI_Aint test_object_ptr;
    MPI_Get_address(test_object, &test_object_ptr);
    MPI_Get_address(&test_object->get_head(), &displacements[0]);
    MPI_Get_address(&test_object->get_tail(), &displacements[1]);
    displacements[0] -= test_object_ptr;
    displacements[1] -= test_object_ptr;
    MPI_Datatype types[2] = {dt_hd.get(), dt_tl.get()};
    MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
    MPI_Type_commit(dt.get_ptr());
  }
  MPI_Datatype get() const {return dt.get();}
};

template <typename Hd>
struct make_mpi_datatype<boost::tuples::cons<Hd, boost::tuples::null_type> > : make_mpi_datatype_base {
  make_mpi_datatype<Hd> dt_hd;
  make_mpi_datatype(): dt_hd() {}
  MPI_Datatype get() const {return dt_hd.get();}
};

#ifndef BOOST_NO_0X_HDR_TUPLE
template <size_t I, size_t N, typename Tuple, typename DTTuple>
struct get_element_addresses_and_types {
  static void go(const Tuple& t, const DTTuple& elts_dt, MPI_Aint addrs[], MPI_Datatype types[]) {
    MPI_Get_address((void*)&std::get<I>(t), &addrs[I]);
    types[I] = std::get<I>(elts_dt).get();
    get_element_addresses_and_types<I + 1, N, Tuple, DTTuple>::go(t, elts_dt, addrs, types);
  }
};

template <size_t N, typename Tuple, typename DTTuple>
struct get_element_addresses_and_types<N, N, Tuple, DTTuple> {
  static void go(const Tuple&, const DTTuple&, MPI_Aint[], MPI_Datatype[]) {}
};

template <typename... Elts>
struct make_mpi_datatype<std::tuple<Elts...>> : make_mpi_datatype_base {
  std::tuple<make_mpi_datatype<Elts> ...> elts_dt;
  scoped_mpi_datatype dt;

  make_mpi_datatype(): elts_dt() {
    int blocklengths[sizeof...(Elts)];
    std::fill(blocklengths, blocklengths + sizeof...(Elts), 1);
    char dummy;
    std::tuple<Elts...> *test_object = (std::tuple<Elts...>*)(&dummy);
    MPI_Aint test_object_ptr;
    MPI_Get_address((void*)test_object, &test_object_ptr);
    MPI_Aint displacements[sizeof...(Elts)];
    MPI_Datatype types[sizeof...(Elts)];
    get_element_addresses_and_types<0, sizeof...(Elts), std::tuple<Elts...>, std::tuple<make_mpi_datatype<Elts>...>>::go(*test_object, elts_dt, displacements, types);
    for (size_t i = 0; i < sizeof...(Elts); ++i) {
      displacements[i] -= test_object_ptr;
    }
    MPI_Type_create_struct(sizeof...(Elts), blocklengths, displacements, types, dt.get_ptr());
    MPI_Type_commit(dt.get_ptr());
  }
  MPI_Datatype get() const {return dt.get();}
};
#endif

extern detail::type_info_map<boost::shared_ptr<make_mpi_datatype_base> > mpi_datatype_map;

template <typename T>
void register_mpi_datatype() {
  mpi_datatype_map.insert(detail::get_type_info<T>(), boost::static_pointer_cast<make_mpi_datatype_base>(boost::make_shared<make_mpi_datatype<T> >()));
}

MPI_Datatype get_mpi_datatype(const std::type_info& ti);

void register_builtin_mpi_datatypes();
void clear_mpi_datatype_map();

}

#endif // AMPLUSPLUS_MAKE_MPI_DATATYPE_HPP
