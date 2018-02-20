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

#ifndef __AGM_HPP
#define __AGM_HPP

#include "../../utils/utils.hpp"
#include "../machine/access_metadata.hpp"

#include <vector>
#include <tuple>

// strict weak ordering
struct chaotic {
public:
  template <typename T>
  bool operator()(T i, T j) {
    return false;
  }
};

struct chaotic_ordering_gen {
public:
  typedef chaotic strict_ordering;
};


template <typename WorkItem,
          typename StrictWeakOrderingRelation,
	  typename MachineModel>
class agm {

  typedef typename MachineModel::buckets_t buckets_t;

private:
  template<typename BucketIterator>
  void shuffle(BucketIterator begin, BucketIterator end) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(begin, end, g);
  }

public:
  agm() {}

  void operator()(MachineModel& machine,
                  std::vector<WorkItem>& initial_workitems) {

    typedef typename buckets_t::bucket_iterator_t workitem_iterator_type;

    debug("Pushing initial workitems ...");
    machine.push_work_items(initial_workitems.begin(), initial_workitems.end());

    while (! machine.all_buckets_empty()) {
      workitem_iterator_type w_begin, w_end;

      // Get an equivalence class of workitems
      typename buckets_t::bucket_t* top_bucket = machine.get_top_bucket();

      while ((top_bucket != NULL) && 
	     (!top_bucket->empty())) {  // for EAGM this would be different

	// get a work-item from the top bucket
	auto w = top_bucket->pop_work_item();

	// send the work-item to appropriate destination
	machine.send(w); 

	// if currently processing bucket is empty, then
	// flush all the queues
	// this is like end of epoch
	if (top_bucket->empty()) {
	  machine.flush_all_queues();
	}
      }

      // We can delete the top bucket
      machine.pop_bucket();

      // We are done processing a bucket.
      // Therefore, synchronize
      machine.synchronize();
    }
  }
};

#endif // __AGM_HPP
