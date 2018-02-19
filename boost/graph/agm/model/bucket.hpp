#ifndef __AGM_UTILS
#define __AGM_UTILS

#include <string>
#include <iostream>
#include <list>
#include <unordered_set>
#include <random>
#include <tuple>
#include <type_traits>
#include <assert.h>
#include <algorithm>

#include <boost/thread/locks.hpp>
#include <boost/parallel/append_buffer.hpp>

#include <am++/make_mpi_datatype.hpp>

#include "bucket_stats.hpp"

namespace boost { namespace graph { namespace agm {

template<typename work_item>
class bucket_structure {
  /**
   * The representation : This work item
   * is used to do the comparison with the strict weak
   * ordering relation. A bucket structure must be created by providing
   * a representation.
   **/

public:

  typedef append_buffer<work_item, 10u> Bucket;

  bucket_structure(work_item& rep) : representation(rep) {
    boost::shared_ptr<Bucket> p(new Bucket);
    buffer.swap(p);
  }

  bool empty() {
    return (buffer->size() == 0);
  }

  void clear() {
    buffer->clear();
  }

  const work_item& get_representation() {
    return representation;
  }

  void push_back(const work_item& wi) {
    buffer->push_back(std::move(wi));
  }

  Bucket* get_bucket() {
    return buffer.get();
  }

private:
  boost::shared_ptr<Bucket> buffer;
  work_item representation;

};

      
template<typename work_item,
         typename strict_weak_ordering>
class buckets {

  // concept assert on whether work_item has < operators defined
public:
  typedef bucket_structure < work_item >  bucket_ds_t;

private:
  std::list< bucket_ds_t* > all_buckets;
  uint64_t index; // only incrementing
  strict_weak_ordering swo;
  int ranks;
  MPI_Datatype dt;
  typename std::list< bucket_ds_t* >::iterator current_bucket;
  bucket_stats stats;
  boost::shared_mutex rwmtx;

  void initialize() {
    current_bucket = all_buckets.begin();
  }
  
public:
  buckets() : index(0),
              ranks(1),
              dt(amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item>())){}

  buckets(strict_weak_ordering& _swo) : index(0),
                                        ranks(1),
                                        dt(amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item>())),
                                        swo(_swo){ initialize(); }

  buckets(int _r,
          strict_weak_ordering& _swo) : index(0),
                                        ranks(_r),
                                        dt(amplusplus::get_mpi_datatype(amplusplus::detail::get_type_info<work_item>())),
                                        swo(_swo){ initialize(); }

  buckets(int _r,
          MPI_Datatype _dt,
          strict_weak_ordering& _swo) : index(0),
                                        ranks(_r),
                                        dt(_dt),
                                        swo(_swo){ initialize(); }
  
  ~buckets() {
    typename std::list< bucket_ds_t* >::iterator ite = all_buckets.begin();
    while (ite != all_buckets.end()) {
      delete (*ite);
      ++ite;
    }

    all_buckets.clear();
  }

  void print_stats() {
    stats.print();
  }

  bool empty() {
    return (all_buckets.size() == 0);
  }

  typename bucket_ds_t::Bucket*
  get_current_bucket() {
    assert(current_bucket != all_buckets.end());
    debug("Inside get_current_bucket ...");
    return (*current_bucket)->get_bucket();
  }

  // No other threads should work while
  // popping the current bucket
  // therefore, we dont need to lock this
  void pop_current_bucket(bool gotobegin = false) {
    assert(all_buckets.size() != 0);
    (*current_bucket)->clear();
    delete (*current_bucket);
    current_bucket = all_buckets.erase(current_bucket);
    // if current_bucket points to end, then point it
    // to the begining
    if (gotobegin ||
        (current_bucket == all_buckets.end()))
      current_bucket = all_buckets.begin();
  }

  // returns true if current bucket is empty in
  // all ranks
  bool synchronize_current_bucket() {
    // is current bucket in local node null ? if not get
    // the representation
    // do a MPI all to all
    //    work_item* recvwis = (work_item*)malloc(sizeof(work_item) * ranks);
    //std::memset((void*)recvwis, 0, sizeof(work_item) * ranks);

    std::vector<work_item> recvwis;
    recvwis.resize(ranks);
    work_item tosend;
    
    if ((current_bucket != all_buckets.end()) &&
        (!(*current_bucket)->empty())) {
      tosend = (*current_bucket)->get_representation();
      // MPI all gather api is confusing !
      // The receive buffer size should indicate the number of
      // elements to be received per rank not the total buffer
      // size !
    } else {
      // Assuming work item first element is a vertex
      std::get<0>(tosend) = UINT64_MAX;
    }

    MPI_Allgather((void*)&tosend, 1, dt, &recvwis[0], 1, dt, MPI_COMM_WORLD);        

    bool all_empty = true;
    for (int i=0; i < ranks; ++i) {
      if (std::get<0>(recvwis[i]) != UINT64_MAX) {
        empty_push(recvwis[i]);
        all_empty = false;
      }     
    }

    return all_empty;
    
  }


  // For debugging
#ifdef PRINT_DEBUG
  void print_wi(pwi w) {
    std::cout << " push v : " << std::get<0>(w) << std::endl;
  }

  void print_wi(int w) {
    std::cout << " push v : " << w << std::endl;
  }
#endif

  void push(work_item& wi, int tid=-1/*not accessed*/) {

    while(true) {
      
      typename std::list< bucket_ds_t* >::iterator ite;
      uint64_t localindex = index;

      // Reads
      {
        boost::shared_lock<boost::shared_mutex> lock(rwmtx);
      
        // first see whether there is a bucket for wi,
        // if so push that to the bucket.
        ite = all_buckets.begin();
        while (ite != all_buckets.end()) {
          bucket_ds_t* pbucket = (*ite);
          // bucket can be empty
          const work_item& representation = pbucket->get_representation();
          // is wi and representation are comparable ?
          if (!swo(wi, representation) && !swo(representation, wi)) {
            // not comparable insert to the current bucket
            pbucket->push_back(wi);
            return;
          } else if (swo(representation, wi)) {
            // if representation < wi, then move to next bucket
            ++ite;
          } else {
            // this means wi < representation, then break
            localindex = index;
            break;
          }
        }
      }

      // Writes
      {
        boost::unique_lock< boost::shared_mutex > lock(rwmtx);        

        if (localindex != index)
          continue;
        
        bool wasempty = false;
        if (all_buckets.size() == 0)
          wasempty = true;
        
        // value change successful
        bucket_ds_t* bucket = new bucket_ds_t(wi);
        bucket->push_back(wi);
        all_buckets.insert(ite, bucket);

        // should be inside an ifdef
        stats.increment_buckets();
        // should be inside an ifdef
        
        if (wasempty)
          current_bucket = all_buckets.begin();

        ++index;
        return;
      }
    }
  }

  //Note : This should only be used while synchronization
  // If we are synchronizing buckets globally
  // then every current_bucket should point to the
  // begining of the bucket list.
  void empty_push(work_item& wi) {

    typename std::list< bucket_ds_t* >::iterator ite;
    // first see whether there is a bucket for wi,
    // if so push that to the bucket.
    ite = current_bucket;
    while (ite != all_buckets.end()) {
      bucket_ds_t* pbucket = (*ite);
      // bucket can be empty
      const work_item& representation = pbucket->get_representation();
      // is wi and representation are comparable ?
      if (!swo(wi, representation) && !swo(representation, wi)) {
        return;
      } else if (swo(representation, wi)) {
        // if representation < wi, then move to next bucket
        ++ite;
      } else {
        // this means wi < representation, then break
        break;
      }
    }

    bool wasempty = false;
    if (all_buckets.size() == 0)
      wasempty = true;
        
    // value change successful
    bucket_ds_t* bucket = new bucket_ds_t(wi);
    all_buckets.insert(ite, bucket);

    // should be inside an ifdef
    stats.increment_buckets();
    // should be inside an ifdef
        
    if (wasempty)
      current_bucket = all_buckets.begin();
  }

  
  template <typename WorkItemIterator>
  void
  push(WorkItemIterator begin, WorkItemIterator end) {
    while(begin != end) {
      push(*begin);
      ++begin;
    }
  }

  template<typename a_bucket_t>
  void print_bucket_elements(a_bucket_t* _bkt) {
    typename a_bucket_t::size_type current_bucket_start, current_bucket_end;
    current_bucket_start = 0;
    current_bucket_end = _bkt->size();

    for (typename a_bucket_t::size_type i = current_bucket_start;
         i < current_bucket_end ; i++) {
      work_item& wi = (*_bkt)[i];
      std::cout << "(" << std::get<0>(wi) << ", " << std::get<1>(wi) << "), ";
    }
  }

  uint64_t get_buckets_created() {
    return stats.get_number_of_bkst_created();
  }

  void print_buckets() {
    stats.print();
    typename std::list< bucket_ds_t* >::iterator ite = all_buckets.begin();
    int j=0;
    for (; ite != all_buckets.end(); ++ite) {
      bucket_ds_t* pbucket = (*ite);
      typename bucket_ds_t::Bucket* theBucket = pbucket->get_bucket();
      std::cout << "Bucket : " << j
                << " [REP : " << "(" << std::get<0>(pbucket->get_representation())
                << ", " << std::get<1>(pbucket->get_representation()) << ")]"
                << ", elements : ";
      print_bucket_elements(theBucket);
      std::cout << std::endl;

      ++j;
    }
  }

};

}}}

#endif
