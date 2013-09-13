#ifndef serial
#include <mpi.h>
#endif

int MPISum_integer(int *global, int *var, int size)
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  return(0);
}

int MPISum_double(double *global, double *var, int size)
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  return(0);
}
