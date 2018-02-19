// Copyright (C) 2014 The Trustees of Indiana University.
 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Andrew Lumsdaine
//           Marcin Zalewski
//           Thejaka Kanewala

#ifndef PBGL2_CALC_TIME
#define PBGL2_CALC_TIME
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

typedef double time_type;

inline time_type get_time()
{
  return MPI_Wtime();
#if 0
  timeval tp;
  gettimeofday(&tp, 0);
  return tp.tv_sec + tp.tv_usec / 1000000.0;
#endif
}
#endif
