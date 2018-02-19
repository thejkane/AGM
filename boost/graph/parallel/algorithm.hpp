// Copyright 2004 and 2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Douglas Gregor
//           Andrew Lumsdaine

#ifndef BOOST_PARALLEL_ALGORITHM_HPP
#define BOOST_PARALLEL_ALGORITHM_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <boost/graph/parallel/graph_utility.hpp>
#include <boost/optional.hpp>
#include <boost/config.hpp> // for BOOST_STATIC_CONSTANT
#include <vector>
#include <functional>

namespace boost { namespace parallel {

  template<typename BinaryOp>
  struct is_commutative
  {
    BOOST_STATIC_CONSTANT(bool, value = false);
  };

  template<typename T>
  struct minimum : std::binary_function<T, T, T>
  {
    const T& operator()(const T& x, const T& y) const { return x < y? x : y; }
  };

  template<typename T>
  struct maximum : std::binary_function<T, T, T>
  {
    const T& operator()(const T& x, const T& y) const { return x < y? y : x; }
  };

  template<typename T>
  struct sum : std::binary_function<T, T, T>
  {
    const T operator()(const T& x, const T& y) const { return x + y; }
  };

  template<typename T, typename BinaryOperation>
  class all_reduce {

    struct reduce_handler {
      reduce_handler() : self(NULL) {}
      reduce_handler(all_reduce& self) : self(&self) {}

      void operator() (amplusplus::rank_type, const T* data, size_t) const { self->handler(data); }

    protected:
      all_reduce* self;
    };

  public:

    all_reduce(amplusplus::transport& trans, BinaryOperation bin_op)
      : transport(trans), bin_op(bin_op), reduce_msg(transport.template create_message_type<T>())
    { 
      reduce_msg.set_max_count(1);
      reduce_msg.set_handler(reduce_handler(*this));
    }
    
    // This can only be called from a single thread, perform the local reduction separately 
    T operator()(const T& in_val) 
    {
      int nthreads = transport.get_nthreads();
      transport.set_nthreads(1);

      out_val = in_val;
      {
	amplusplus::scoped_epoch epoch(transport);
	
	// TODO: Use reductions to perform the bin_op's in the network
	for (typename amplusplus::transport::rank_type i = (transport.rank() + 1) % transport.size(); 
	     i != transport.rank(); i = (i + 1) % transport.size()) {
	  reduce_msg.message_being_built(i);
	  reduce_msg.send(&in_val, 1, i, empty_deleter());
	}
      }

      transport.set_nthreads(nthreads);

      return out_val;
    }

    void handler(const T* data) { out_val = bin_op(out_val, *data); }

  private:
    
    amplusplus::transport& transport;
    T out_val;
    BinaryOperation bin_op;
    amplusplus::message_type<T> reduce_msg;

    struct empty_deleter { void operator()() const {} };
  };

  // TODO: Implement other collectives... scan, etc.

  // Sum over unsigned long long using scoped epoch value
  long unsigned int all_reduce_sum(amplusplus::transport& transport, const long unsigned int& val) {

    long unsigned int result;

    { amplusplus::scoped_epoch_value epoch(transport, val, result); }

    return result;
  }

  // all_gather
  template <typename T, typename MessageGenerator = 
                          amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
  class all_gather {

    // For explicit addressing
    typedef std::pair<amplusplus::transport::rank_type, T> msg_data_type;

    struct gather_handler {
      gather_handler() : self(NULL) {}
      gather_handler(all_gather& self) : self(&self) {}

      void operator() (const msg_data_type& data) const { self->handler(data.second); }

    protected:
      all_gather* self;
    };

    typedef typename MessageGenerator::template call_result<msg_data_type, gather_handler, map_of_1st,
				       const amplusplus::duplicate_removal_t<boost::parallel::identity<T> > >::type 
      GatherMessage;

  public:

    all_gather(amplusplus::transport& trans,
               MessageGenerator message_gen = 
                 MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))

      : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<msg_data_type>(), 0)),
        transport(trans),
 	gather_msg(message_gen, transport, map_of_1st(), amplusplus::duplicate_removal(boost::parallel::identity<T>()))
    { 
      gather_msg.set_handler(gather_handler(*this));
    }

    // This can only be called from a single thread, perform the local reduction separately 
    void operator()(std::vector<T>& in, std::vector<T>& out)
    {
      int nthreads = transport.get_nthreads();
      transport.set_nthreads(1);

      out_ptr = &out;
      out.insert(out.end(), in.begin(), in.end());
      {
	amplusplus::scoped_epoch epoch(transport);
	
	// TODO: Use reductions to perform the bin_op's in the network
	for (typename amplusplus::transport::rank_type i = (transport.rank() + 1) % transport.size(); 
	     i != transport.rank(); i = (i + 1) % transport.size())
	  for (size_t j = 0 ; j < in.size() ; ++j)  
	    gather_msg.send(std::make_pair(i, in[j])); // NGE: Should send whole in vector rather than single elements here
   	    // TODO: Multicast would be nice here!
	    // NGE: Explicitly including the destination for every element here is awful
	    //      we really need an explicit send-to-rank in OBA 
      }

      transport.set_nthreads(nthreads);

      return;
    }

    void handler(const T& data) {
      out_ptr->insert(out_ptr->end(), data);
    }

  private:

    const int dummy_first_member_for_init_order; // Unused    

    amplusplus::transport& transport;
    std::vector<T>* out_ptr;
    GatherMessage gather_msg;
  };

} } // end namespace boost::parallel

#endif // BOOST_PARALLEL_ALGORITHM_HPP
