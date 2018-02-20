#ifndef NBC_STUB_H
#define NBC_STUB_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#define NBC_Request MPI_Request

#define NBC_Ibarrier MPI_Ibarrier
#define NBC_Iallreduce MPI_Iallreduce
#define NBC_Wait MPI_Wait

static inline int NBC_Test(NBC_Request* req, int* flag, MPI_Status* status) {
  int res = MPI_Test(req, flag, status);
  if (res == MPI_SUCCESS) *flag = (*flag ? 0 : 3); // Change true to NBC_OK and false to NBC_CONTINUE
  return res;
}

#define NBC_OK MPI_SUCCESS

#ifdef __cplusplus
}
#endif

#endif /* NBC_STUB_H */
