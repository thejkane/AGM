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

#ifndef AMPLUSPLUS_HPP
#define AMPLUSPLUS_HPP

#include <stdint.h>
#include <vector>
#include <boost/utility.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread.hpp>
#include <set>
#include <functional>
#include <am++/traits.hpp>
#include <am++/transport.hpp>
#include <am++/termination_detector.hpp>
#include <am++/scoped_epoch.hpp>
#include <am++/message_type_generators.hpp>

// Create AM transport in one thread
// Begin/end epochs in all threads that will be doing progress
//   Term. detection is only done once all threads have ended the epoch, but
//   they all run progress until tdetect actually starts
// Sends and non-epoch-ending progress calls can happen from anywhere
//   Progress must be called from multiple threads to parallelize handlers (or
//   the user code/coalescing layer must dispatch handler calls itself)

// Note: termination detectors are multi-threaded, including reset() and
// finish()

// Fast coalesced send:
// Outgoing message buffer is two contiguous words (for DWCAS)
//   Message count (or # remaining)
//   Buffer ptr (or index into table)
// Send does DWCAS to increment msg count, or if the buffer is going to be
// full, swap in another buffer and zero as the count
// Another counter in each buffer for number of entries fully written
// Send is:
//   Atomically read and update (count, ptr) to
//     (count + 1 == size) ? (0, newptr) : (count + 1, ptr)
//   If first branch was taken, wait until ptr->written_count == size then send
//   buffer, replacing newptr with a new buffer from the pool
//   If second branch was taken, write element count of ptr->data and increment
//   ptr->written_count
// Each thread has a new buffer in TLS for fast access for sends w/o hitting
// the buffer pool

#endif // AMPLUSPLUS_HPP
