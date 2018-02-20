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

#ifndef AMPLUSPLUS_SCOPED_EPOCH_HPP
#define AMPLUSPLUS_SCOPED_EPOCH_HPP

#include <boost/utility.hpp>
#include <boost/noncopyable.hpp>
#include <am++/transport.hpp>

namespace amplusplus {

class scoped_epoch: boost::noncopyable {
  transport tr;
  public:
  scoped_epoch(transport tr): tr(tr) {tr.begin_epoch();}
  ~scoped_epoch() {tr.end_epoch();}
};

class scoped_epoch_value: boost::noncopyable {
  transport tr;
  const unsigned long& read_value;
  unsigned long& sum;
  public:
  scoped_epoch_value(transport tr,
                     const unsigned long& read_value, unsigned long& sum)
    : tr(tr), read_value(read_value), sum(sum) {tr.begin_epoch();}
  ~scoped_epoch_value() {sum = tr.end_epoch_with_value(read_value);}
};

}

#endif // AMPLUSPLUS_SCOPED_EPOCH_HPP
