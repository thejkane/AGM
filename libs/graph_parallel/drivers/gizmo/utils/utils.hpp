// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

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

// The final execution time in terms of cycles 
// This is a global variable
extern uint64_t execution_time;

/*template<typename CBClass>
class Callback {

  typedef Callback<CBClass> self_type;

private:
  CBClass& cb;

public:
  Callback(CBClass& _cb) : cb(_cb) {}

  void static Call(void* po, uint64_t t) {
    self_type* st = (self_type*)po;
    (*st)(t);
  }

  template<typename arg_t>
  void operator()(arg_t arg) {
    cb.flush(arg);
  }
  };*/


void invoke_gizmo_cb(void* _m, uint64_t t);

// class encapsulate some utility functions
class util {

private:
  static void* p_callback;

  static void time(uint64_t val, bool force_flush=false) {
    if (val > execution_time) {
      execution_time = val; 
      assert(p_callback != NULL);
      invoke_gizmo_cb(p_callback, execution_time);
    } else {
      if (force_flush) {
	assert(p_callback != NULL);
	invoke_gizmo_cb(p_callback, execution_time);
      }
    }
  }

public:
  
  void static set_callback(void* _pcb) {
    assert(_pcb != NULL);
    p_callback = _pcb;
  }

  void static invoke_call_back() {
    invoke_gizmo_cb(p_callback, execution_time);
  }
  
  template<typename work_item_t>
  void static increment_time(work_item_t& wi) {
    auto val = ++(std::get<2>(wi));     
    time(val);
  }

  template<typename work_item_t>
  void static increment_time_no_cb(work_item_t& wi) {
    uint64_t val = ++(std::get<2>(wi));     
    if (val > execution_time) {
      execution_time = val; 
    }
  }

  template<typename work_item_t>
  void static increment_time(work_item_t& wi, uint64_t t) {
    std::get<2>(wi) = std::get<2>(wi) + t;
    auto val = std::get<2>(wi);
    time(val);
  }

  void static increment_time(uint64_t val, bool force_flush=false) {
    time(val, force_flush);
  }
};

void info(const char* str) {
  //std::cout << "[INFO] " << str << std::endl;
}

void error(const char* str) {
  std::cout << "[ERROR] " << str << std::endl;
}

void debug(const char* str) {
  //std::cout << "[DEBUG] " << str << std::endl;
}

void debug(int str) {
  //std::cout << "[DEBUG] " << str << std::endl;
}

//mt19937 gen( chrono::system_clock::now().time_since_epoch().count() );

template<typename work_item>
class BucketStructure {
  /**
   * The representation : This work item
   * is used to do the comparison with the strict weak
   * ordering relation. A bucket structure must be created by providing
   * a representation.
   **/

public:

  typedef typename std::vector < work_item >::iterator iterator;

  BucketStructure(work_item& rep) : representation(rep) {}

  bool empty() {
    return bucket.empty();
  }

  void shuffle() {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(bucket.begin(), bucket.end(), g);
  }

  iterator begin() {
    return bucket.begin();
  }

  iterator end() {
    return bucket.end();
  }

  void clear() {
    bucket.clear();
  }

  work_item get_representation() {
    return representation;
  }

  void push_back(const work_item& wi) {
    bucket.push_back(wi);
  }

  work_item pop_work_item() {
    // caller must first check whether bucket is empty
    // if bucket is not empty caller should invoke this function
    assert(!bucket.empty());
    shuffle();

    auto it = bucket.begin();
    work_item w = (*it);
    bucket.erase(it);
    return w;
  }

private:
  std::vector < work_item >  bucket;
  work_item representation;

};

template<typename work_item,
         typename strict_weak_ordering>
class buckets {

  // concept assert on whether work_item has < operators defined
public:
  typedef BucketStructure < work_item >  bucket_t;

private:
  std::list< bucket_t* > all_buckets;
  strict_weak_ordering swo;

public:
  typedef typename bucket_t::iterator bucket_iterator_t;

  ~buckets() {
    typename std::list< bucket_t* >::iterator ite = all_buckets.begin();
    while (ite != all_buckets.end()) {
      delete (*ite);
      ++ite;
    }

    all_buckets.clear();
  }

  bool empty() {
    return (all_buckets.size() == 0);
  }

  bucket_t* get_top_bucket() {
    return all_buckets.front();
  }

  std::pair<typename bucket_t::iterator, typename bucket_t::iterator>
  top_bucket() {
    bucket_t* pbucket = get_top_bucket();
    
    // top bucket could be empty
    //assert(!pbucket->empty());
    pbucket->shuffle();

    return std::make_pair(pbucket->begin(), pbucket->end());
  }

  void pop_bucket() {
    assert(all_buckets.size() != 0);
    typename std::list< bucket_t* >::iterator ite = all_buckets.begin();
    (*ite)->clear();
    delete (*ite);
    all_buckets.erase(ite);
  }

  template<typename pwi>
  void print_wi(pwi w) {
    std::cout << " push v : " << std::get<0>(w) << std::endl;
  }

  void print_wi(int w) {
    std::cout << " push v : " << w << std::endl;
  }

  void push(work_item& wi) {

#ifdef PRINT_DEBUG
    print_wi(wi);
#endif

    // if buckets are empty, create new one
    if (all_buckets.size() == 0) {
      // when creating the bucket, pass the current work item
      // as the representing work item
      bucket_t* bucket = new bucket_t(wi);
      bucket->push_back(wi);
      all_buckets.push_back(bucket);
      return;
    }

    typename std::list< bucket_t* >::iterator ite = all_buckets.begin();
    while (ite != all_buckets.end()) {
      bucket_t* pbucket = (*ite);
      // bucket can be empty
      work_item representation = pbucket->get_representation();
      // is wi and representation are comparable ?
      if (!swo(wi, representation) && !swo(representation, wi)) {
        // not comparable insert to the current bucket
        pbucket->push_back(wi);
        return;
      }

      // work items are comparable
      if (swo(wi, representation)) {
        // we need to insert a new bucket before current bucket
        bucket_t* bucket = new bucket_t(wi);
        bucket->push_back(wi);
        all_buckets.insert(ite, bucket);
        return;
      }
      // if representation < wi, then move to next bucket
      ++ite;
    }

    // not yet inserted -- that means no bucket is matching
    // we need to insert a new bucket
    bucket_t* bucket = new bucket_t(wi);
    bucket->push_back(wi);
    all_buckets.push_back(bucket);
  }

  template <typename WorkItemIterator>
  void
  push(WorkItemIterator begin, WorkItemIterator end) {
    while(begin != end) {
      push(*begin);
      ++begin;
    }
  }


  template <typename WorkItemIterator>
  typename bucket_t::iterator
  push(WorkItemIterator begin, WorkItemIterator end,
       typename bucket_t::iterator current) {

    push(begin, end);

    typename std::list< bucket_t* >::iterator ite = all_buckets.begin();
    bucket_t* pbucket = (*ite);
    assert(!pbucket->empty());

    return pbucket->end();
  }


  void print_buckets() {
    typename std::list< bucket_t* >::iterator ite = all_buckets.begin();
    int j=0;
    for (; ite != all_buckets.end(); ++ite) {
      bucket_t* pbucket = (*ite);
      typename bucket_t::iterator itebkt = pbucket->begin();
      for (; itebkt != pbucket->end(); ++itebkt) {
        std::cout << "Bucket : " << j << ", element : " << (*itebkt) << std::endl;
      }
      ++j;
    }
  }

};


#endif
