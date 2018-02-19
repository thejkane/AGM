// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Peter Gottschling
//           Andrew Lumsdaine
#ifndef BOOST_PARALLEL_DISTRIBUTION_HPP
#define BOOST_PARALLEL_DISTRIBUTION_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <cstddef>
#include <vector>
#include <algorithm>
#include <numeric>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/shared_ptr.hpp>
#include <typeinfo>

namespace boost { namespace parallel {

template<typename SizeType = std::size_t>
class variant_distribution
{
public:
  typedef typename amplusplus::transport::rank_type rank_type;
  typedef std::size_t process_size_type;
  typedef SizeType size_type;

private:
  struct basic_distribution
  {
    virtual ~basic_distribution() {}
    virtual size_type block_size(rank_type, size_type) const = 0;
    virtual rank_type in_process(size_type) const = 0;
    virtual size_type local(size_type) const = 0;
    virtual size_type global(size_type) const = 0;
    virtual size_type global(rank_type, size_type) const = 0;
    virtual void* address() = 0;
    virtual const void* address() const = 0;
    virtual const std::type_info& type() const = 0;
  };

  template<typename Distribution>
  struct poly_distribution : public basic_distribution
  {
    explicit poly_distribution(const Distribution& distribution)
      : distribution_(distribution) { }

    virtual size_type block_size(rank_type id, size_type n) const
    { return distribution_.block_size(id, n); }

    virtual rank_type in_process(size_type i) const
    { return distribution_(i); }

    virtual size_type local(size_type i) const
    { return distribution_.local(i); }

    virtual size_type global(size_type n) const
    { return distribution_.global(n); }

    virtual size_type global(rank_type id, size_type n) const
    { return distribution_.global(id, n); }

    virtual void* address() { return &distribution_; }
    virtual const void* address() const { return &distribution_; }
    virtual const std::type_info& type() const { return typeid(Distribution); }

  private:
    Distribution distribution_;
  };

public:
  variant_distribution() { }

  template<typename Distribution>
  variant_distribution(const Distribution& distribution)
    : distribution_(new poly_distribution<Distribution>(distribution)) { }

  size_type block_size(rank_type id, size_type n) const
  { return distribution_->block_size(id, n); }
  
  rank_type operator()(size_type i) const
  { return distribution_->in_process(i); }
  
  size_type local(size_type i) const
  { return distribution_->local(i); }
  
  size_type global(size_type n) const
  { return distribution_->global(n); }

  size_type global(rank_type id, size_type n) const
  { return distribution_->global(id, n); }

  operator bool() const { return distribution_; }

  void clear() { distribution_.reset(); }

  template<typename T>
  T* as()
  {
    if (distribution_->type() == typeid(T))
      return static_cast<T*>(distribution_->address());
    else
      return 0;
  }

  template<typename T>
  const T* as() const
  {
    if (distribution_->type() == typeid(T))
      return static_cast<T*>(distribution_->address());
    else
      return 0;
  }

private:
  shared_ptr<basic_distribution> distribution_;
};

template<typename SizeType = std::size_t>
struct block
{
  explicit block(const amplusplus::transport& trans, SizeType n) 
    : id(trans.rank()), p(trans.size()), n(n) { }

  // If there are n elements in the distributed data structure, returns the number of elements stored locally.
  SizeType block_size(SizeType n) const
  { return (n / p) + ((SizeType)(n % p) > id? 1 : 0); }

  // If there are n elements in the distributed data structure, returns the number of elements stored on processor ID
  template<typename ProcessID>
  SizeType block_size(ProcessID id, SizeType n) const
  { return (n / p) + ((ProcessID)(n % p) > id? 1 : 0); }

  // Returns the processor on which element with global index i is stored
  SizeType operator()(SizeType i) const
  { 
    SizeType cutoff_processor = n % p;
    SizeType cutoff = cutoff_processor * (n / p + 1);

    if (i < cutoff) return i / (n / p + 1);
    else return cutoff_processor + (i - cutoff) / (n / p);
  }

  // Find the starting index for processor with the given id
  template<typename ID>
  SizeType start(ID id) const
  {
    SizeType estimate = id * (n / p + 1);
    ID cutoff_processor = n % p;
    if (id < cutoff_processor) return estimate;
    else return estimate - (id - cutoff_processor);
  }

  // Find the local index for the ith global element
  SizeType local(SizeType i) const
  { 
    SizeType owner = (*this)(i);
    return i - start(owner);
  }

  // Returns the global index of local element i
  SizeType global(SizeType i) const
  { return global(id, i); }

  // Returns the global index of the ith local element on processor id
  template<typename ProcessID>
  SizeType global(ProcessID id, SizeType i) const
  { return i + start(id); }

 private:
  SizeType id; //< The ID number of this processor
  SizeType p;  //< The number of processors
  SizeType n;  //< The size of the problem space
};

// Block distribution with arbitrary block sizes

template<typename SizeType = std::size_t>
struct uneven_block
{
  typedef typename std::vector<SizeType>    size_vector;

  explicit uneven_block(const amplusplus::transport& trans, const std::vector<SizeType>& local_sizes) 
    : id(trans.rank()), p(trans.size()), local_sizes(local_sizes)
  { 
    assert(local_sizes.size() == p);
    local_starts.resize(p + 1);
    local_starts[0] = 0;
    std::partial_sum(local_sizes.begin(), local_sizes.end(), &local_starts[1]);
    n = local_starts[p];
  }

  // To do maybe: enter local size in each process and gather in constructor (much handier)
  // explicit uneven_block(const amplusplus::transport& trans, std::size_t my_local_size) 

  // If there are n elements in the distributed data structure, returns the number of elements stored locally.
  SizeType block_size(SizeType) const
  { return local_sizes[id]; }

  // If there are n elements in the distributed data structure, returns the number of elements stored on processor ID
  template<typename ProcessID>
  SizeType block_size(ProcessID id, SizeType) const
  { return local_sizes[id]; }

  // Returns the processor on which element with global index i is stored
  SizeType operator()(SizeType i) const
  {
    assert (i >= (SizeType) 0 && i < (SizeType) n); // check for valid range
    typename size_vector::const_iterator lb 
      = std::lower_bound(local_starts.begin(), local_starts.end(), (SizeType) i);
    return ((SizeType)(*lb) == i ? lb : --lb) - local_starts.begin();
  }

  // Find the starting index for processor with the given id
  template<typename ID>
  SizeType start(ID id) const 
  {
    return local_starts[id];
  }

  // Find the local index for the ith global element
  SizeType local(SizeType i) const
  { 
    SizeType owner = (*this)(i);
    return i - start(owner);
  }

  // Returns the global index of local element i
  SizeType global(SizeType i) const
  { return global(id, i); }

  // Returns the global index of the ith local element on processor id
  template<typename ProcessID>
  SizeType global(ProcessID id, SizeType i) const
  { return i + start(id); }

 private:
  SizeType              id;           //< The ID number of this processor
  SizeType              p;            //< The number of processors
  SizeType              n;            //< The size of the problem space
  std::vector<SizeType> local_sizes;  //< The sizes of all blocks
  std::vector<SizeType> local_starts; //< Lowest global index of each block
};

template<typename SizeType = std::size_t>
struct oned_block_cyclic
{
  explicit oned_block_cyclic(const amplusplus::transport& trans, SizeType size)
    : id(trans.rank()), p(trans.size()), size(size) { }
      
  SizeType block_size(SizeType n) const
  { 
    return block_size(id, n);
  }

  template<typename ProcessID>
  SizeType block_size(ProcessID i, SizeType n) const
  {
    //    std::cout << "block_size - ProcessID : " << i << " n : " << n << std::endl;
    SizeType all_blocks = n / size;
    SizeType extra_elements = n % size;
    SizeType everyone_gets = all_blocks / p;
    SizeType extra_blocks = all_blocks % p;
    SizeType my_blocks = everyone_gets + (i < extra_blocks? 1 : 0);
    SizeType my_elements = my_blocks * size 
                         + (i == extra_blocks? extra_elements : 0);

    //    std::cout << "my_elements : " << my_elements << std::endl;
    return my_elements;
  }

  SizeType operator()(SizeType i) const
  { 
    auto val = (i / size) % p;
    //std::cout << "operator()" << "i : " << i << " size : " << size << " p : " << p << " val : " << val << std::endl; 
    return val;
  }

  SizeType local(SizeType i) const
  { 
    return ((i / size) / p) * size + i % size;
  }

  SizeType global(SizeType i) const
  { return global(id, i); }

  template<typename ProcessID>
  SizeType global(ProcessID id, SizeType i) const
  { 
    return ((i / size) * p + id) * size + i % size;
  }

 private:
  SizeType id;                   //< The ID number of this processor
  SizeType p;                    //< The number of processors
  SizeType size;                 //< Block size
};

template<typename SizeType = std::size_t>
struct twod_block_cyclic
{
  explicit twod_block_cyclic(const amplusplus::transport& trans,
                             SizeType block_rows, SizeType block_columns,
                             SizeType data_columns_per_row)
    : id(trans.rank()), p(trans.size()), 
      block_rows(block_rows), block_columns(block_columns), 
      data_columns_per_row(data_columns_per_row)
  { }
      
  SizeType block_size(SizeType n) const
  { 
    return block_size(id, n);
  }

  template<typename ProcessID>
  SizeType block_size(ProcessID id, SizeType n) const
  {
    // TBD: This is really lame :)
    int result = -1;
    while (n > 0) {
      --n;
      if ((*this)(n) == id && (int)local(n) > result) result = local(n);
    }
    ++result;

    //    std::cerr << "Block size of id " << id << " is " << result << std::endl;
    return result;
  }

  SizeType operator()(SizeType i) const
  { 
    SizeType result = get_block_num(i) % p;
    //    std::cerr << "Item " << i << " goes on processor " << result << std::endl;
    return result;
  }

  SizeType local(SizeType i) const
  { 
    // Compute the start of the block
    SizeType block_num = get_block_num(i);
    //    std::cerr << "Item " << i << " is in block #" << block_num << std::endl;

    SizeType local_block_num = block_num / p;
    SizeType block_start = local_block_num * block_rows * block_columns;

    // Compute the offset into the block 
    SizeType data_row = i / data_columns_per_row;
    SizeType data_col = i % data_columns_per_row;
    SizeType block_offset = (data_row % block_rows) * block_columns 
                             + (data_col % block_columns);    

    //    std::cerr << "Item " << i << " maps to local index " << block_start+block_offset << std::endl;
    return block_start + block_offset;
  }

  SizeType global(SizeType i) const
  { 
    // Compute the (global) block in which this element resides
    SizeType local_block_num = i / (block_rows * block_columns);
    SizeType block_offset = i % (block_rows * block_columns);
    SizeType block_num = local_block_num * p + id;

    // Compute the position of the start of the block (globally)
    SizeType block_start = block_num * block_rows * block_columns;

//     std::cerr << "Block " << block_num << " starts at index " << block_start
//               << std::endl;

    // Compute the row and column of this block
    SizeType block_row = block_num / (data_columns_per_row / block_columns);
    SizeType block_col = block_num % (data_columns_per_row / block_columns);

    SizeType row_in_block = block_offset / block_columns;
    SizeType col_in_block = block_offset % block_columns;

//     std::cerr << "Local index " << i << " is in block at row " << block_row
//               << ", column " << block_col << ", in-block row " << row_in_block
//               << ", in-block col " << col_in_block << std::endl;

    SizeType result = block_row * block_rows + block_col * block_columns
                    + row_in_block * block_rows + col_in_block;

//     std::cerr << "global(" << i << "@" << id << ") = " << result 
//               << " =? " << local(result) << std::endl;
    assert(i == local(result));
    return result;
  }

 private:

  SizeType get_block_num(SizeType i) const
  {
    SizeType data_row = i / data_columns_per_row;
    SizeType data_col = i % data_columns_per_row;
    SizeType block_row = data_row / block_rows;
    SizeType block_col = data_col / block_columns;
    SizeType blocks_in_row = data_columns_per_row / block_columns;
    SizeType block_num = block_col * blocks_in_row + block_row;
    return block_num;
  }

  SizeType id;                   //< The ID number of this processor
  SizeType p;                    //< The number of processors
  SizeType block_rows;           //< The # of rows in each block
  SizeType block_columns;        //< The # of columns in each block
  SizeType data_columns_per_row; //< The # of columns per row of data
};

template<typename SizeType = std::size_t>
class twod_random
{
  template<typename RandomNumberGen>
  struct random_int
  {
    explicit random_int(RandomNumberGen& gen) : gen(gen) { }

    template<typename T>
    T operator()(T n) const
    {
      uniform_int<T> distrib(0, n-1);
      return distrib(gen);
    }

  private:
    RandomNumberGen& gen;
  };
  
 public:
  template<typename RandomNumberGen>
  explicit twod_random(const amplusplus::transport& trans,
                       SizeType block_rows, SizeType block_columns,
                       SizeType data_columns_per_row,
                       SizeType n,
                       RandomNumberGen& gen)
    : id(trans.rank()), p(trans.size()), 
      block_rows(block_rows), block_columns(block_columns), 
      data_columns_per_row(data_columns_per_row),
      global_to_local(n / (block_rows * block_columns))
  { 
    std::copy(make_counting_iterator(SizeType(0)),
              make_counting_iterator(global_to_local.size()),
              global_to_local.begin());

    random_int<RandomNumberGen> rand(gen);
    std::random_shuffle(global_to_local.begin(), global_to_local.end(), rand);
  }
      
  SizeType block_size(SizeType n) const
  { 
    return block_size(id, n);
  }

  template<typename ProcessID>
  SizeType block_size(ProcessID id, SizeType n) const
  {
    // TBD: This is really lame :)
    int result = -1;
    while (n > 0) {
      --n;
      if ((*this)(n) == id && (int)local(n) > result) result = local(n);
    }
    ++result;

    //    std::cerr << "Block size of id " << id << " is " << result << std::endl;
    return result;
  }

  SizeType operator()(SizeType i) const
  { 
    SizeType result = get_block_num(i) % p;
    //    std::cerr << "Item " << i << " goes on processor " << result << std::endl;
    return result;
  }

  SizeType local(SizeType i) const
  { 
    // Compute the start of the block
    SizeType block_num = get_block_num(i);
    //    std::cerr << "Item " << i << " is in block #" << block_num << std::endl;

    SizeType local_block_num = block_num / p;
    SizeType block_start = local_block_num * block_rows * block_columns;

    // Compute the offset into the block 
    SizeType data_row = i / data_columns_per_row;
    SizeType data_col = i % data_columns_per_row;
    SizeType block_offset = (data_row % block_rows) * block_columns 
                             + (data_col % block_columns);    

    //    std::cerr << "Item " << i << " maps to local index " << block_start+block_offset << std::endl;
    return block_start + block_offset;
  }

 private:

  SizeType get_block_num(SizeType i) const
  {
    SizeType data_row = i / data_columns_per_row;
    SizeType data_col = i % data_columns_per_row;
    SizeType block_row = data_row / block_rows;
    SizeType block_col = data_col / block_columns;
    SizeType blocks_in_row = data_columns_per_row / block_columns;
    SizeType block_num = block_col * blocks_in_row + block_row;
    return global_to_local[block_num];
  }

  SizeType id;                   //< The ID number of this processor
  SizeType p;                    //< The number of processors
  SizeType block_rows;           //< The # of rows in each block
  SizeType block_columns;        //< The # of columns in each block
  SizeType data_columns_per_row; //< The # of columns per row of data
  std::vector<SizeType> global_to_local;
};

template<typename SizeType = std::size_t>
class random_distribution
{
  template<typename RandomNumberGen>
  struct random_int
  {
    explicit random_int(RandomNumberGen& gen) : gen(gen) { }

    template<typename T>
    T operator()(T n) const
    {
      uniform_int<T> distrib(0, n-1);
      return distrib(gen);
    }

  private:
    RandomNumberGen& gen;
  };
  
 public:
  template<typename RandomNumberGen>
  random_distribution(const amplusplus::transport& trans, RandomNumberGen& gen,
                      SizeType n)
    : base(trans, n), local_to_global(n), global_to_local(n)
  {
    std::copy(make_counting_iterator(SizeType(0)),
              make_counting_iterator(n),
              local_to_global.begin());

    random_int<RandomNumberGen> rand(gen);
    std::random_shuffle(local_to_global.begin(), local_to_global.end(), rand);
                        

    for (typename std::vector<SizeType>::size_type i = 0; i < n; ++i)
      global_to_local[local_to_global[i]] = i;
  }

  SizeType block_size(SizeType n) const
  { return base.block_size(n); }

  template<typename ProcessID>
  SizeType block_size(ProcessID id, SizeType n) const
  { return base.block_size(id, n); }

  SizeType operator()(SizeType i) const
  {
    return base(global_to_local[i]);
  }

  SizeType local(SizeType i) const
  { 
    return base.local(global_to_local[i]);
  }

  template<typename ProcessID>
  SizeType global(ProcessID p, SizeType i) const
  { 
    return local_to_global[base.global(p, i)];
  }

  SizeType global(SizeType i) const
  { 
    return local_to_global[base.global(i)];
  }

 private:
  block<SizeType> base;
  std::vector<SizeType> local_to_global;
  std::vector<SizeType> global_to_local;
};

} } // end namespace boost::parallel

#endif // BOOST_PARALLEL_DISTRIBUTION_HPP


//  LocalWords:  SizeType
