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

#ifndef AMPLUSPLUS_DETAIL_ID_ASSIGNER_HPP
#define AMPLUSPLUS_DETAIL_ID_ASSIGNER_HPP

namespace amplusplus {
  namespace detail {

class id_assigner {
  unsigned int highest_ever_assigned_plus_one;
  std::vector<unsigned int> free_values;

  public:
  id_assigner(): highest_ever_assigned_plus_one(0) {}
  unsigned int allocate() {
    if (!free_values.empty()) {
      unsigned int v = free_values.back();
      free_values.pop_back();
      return v;
    }
    return highest_ever_assigned_plus_one++;
  }
  void free(unsigned int v) {
    if (v + 1 == highest_ever_assigned_plus_one) {
      --highest_ever_assigned_plus_one;
    } else {
      free_values.push_back(v);
    }
  }
};

class scoped_id {
  id_assigner& assigner;
  unsigned int value;

  public:
  scoped_id(id_assigner& assigner): assigner(assigner) {
    value = assigner.allocate();
  }
  ~scoped_id() {
    assigner.free(value);
  }
  unsigned int get_value() const {return value;}
};

  }
}

#endif // AMPLUSPLUS_DETAIL_ID_ASSIGNER_HPP
