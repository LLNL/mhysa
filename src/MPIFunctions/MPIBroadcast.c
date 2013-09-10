#include <mpivars.h>

int MPIBroadcast_double(double *x, int size, int root)
{
#ifndef serial
  MPI_Bcast(x,size,MPI_DOUBLE,root,MPI_COMM_WORLD);
#endif
  return(0);
}

int MPIBroadcast_integer(int *x, int size, int root)
{
#ifndef serial
  MPI_Bcast(x,size,MPI_INT,root,MPI_COMM_WORLD);
#endif
  return(0);
}

int MPIBroadcast_character(char *x, int size, int root)
{
#ifndef serial
  MPI_Bcast(x,size,MPI_CHAR,root,MPI_COMM_WORLD);
#endif
  return(0);
}
