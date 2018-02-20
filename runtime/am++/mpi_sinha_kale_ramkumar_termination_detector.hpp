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

#ifndef AMPLUSPLUS_MPI_SINHA_KALE_RAMKUMAR_TERMINATION_DETECTOR_HPP
#define AMPLUSPLUS_MPI_SINHA_KALE_RAMKUMAR_TERMINATION_DETECTOR_HPP

#include <am++/transport.hpp>
#include <am++/termination_detector.hpp>

namespace amplusplus {
  // Unbounded depth
  // From http://charm.cs.uiuc.edu/papers/QuiescenceINTL94.pdf

  // Note: initialize() is not thread-safe; everything else is
  termination_detector make_mpi_sinha_kale_ramkumar_termination_detector(transport& trans);
}

#endif // AMPLUSPLUS_MPI_SINHA_KALE_RAMKUMAR_TERMINATION_DETECTOR_HPP
