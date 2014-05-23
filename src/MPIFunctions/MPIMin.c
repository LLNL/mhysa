#include <mpivars.h>

int MPIMin_integer(int *global, int *var, int size,void *comm)
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_INT,MPI_MIN,*((MPI_Comm*)comm));
#endif
  return(0);
}

int MPIMin_double(double *global, double *var, int size,void *comm)
{
#ifdef serial
  int i;
  for (i = 0; i < size; i++)  global[i] = var[i];
#else
  MPI_Allreduce((var==global?MPI_IN_PLACE:var),global,size,MPI_DOUBLE,MPI_MIN,*((MPI_Comm*)comm));
#endif
  return(0);
}
