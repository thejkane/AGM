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

// Strange file name is so Makefile system recognizes that this file needs to
// be built with mpic++.

#include <config.h>

#include <am++/make_mpi_datatype.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <boost/smart_ptr.hpp>

namespace amplusplus {
  detail::type_info_map<boost::shared_ptr<make_mpi_datatype_base> > mpi_datatype_map;

MPI_Datatype get_mpi_datatype(const std::type_info& ti) {
  const boost::shared_ptr<make_mpi_datatype_base>* p = mpi_datatype_map.lookup(ti);
  if (!p) {
    fprintf(stderr, "Did not find MPI datatype for C++ type %s; make sure to register it with amplusplus::register_mpi_datatype()\n",
            ti.name());
    abort();
  }
  return (*p)->get();
}

void register_builtin_mpi_datatypes() {
  register_mpi_datatype<char>();
  register_mpi_datatype<signed char>();
  register_mpi_datatype<unsigned char>();
  register_mpi_datatype<short>();
  register_mpi_datatype<unsigned short>();
  register_mpi_datatype<int>();
  register_mpi_datatype<unsigned int>();
  register_mpi_datatype<long>();
  register_mpi_datatype<unsigned long>();
  register_mpi_datatype<long long>();
  register_mpi_datatype<unsigned long long>();
  register_mpi_datatype<float>();
  register_mpi_datatype<double>();
  register_mpi_datatype<long double>();
}

void clear_mpi_datatype_map() {
  mpi_datatype_map.clear();
}

}
