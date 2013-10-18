#include <stdlib.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPICreateCommunicators(int ndims,void *m)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          i,n,color,key;
  int          *ip,*iproc;

  mpi->comm = (MPI_Comm*) calloc (ndims, sizeof(MPI_Comm));
  if (ndims == 1) MPI_Comm_dup(mpi->world,mpi->comm);
  else {
    ip    = (int*) calloc (ndims-1,sizeof(int));
    iproc = (int*) calloc (ndims-1,sizeof(int));
    for (n=0; n<ndims; n++) {
      int tick=0; 
      for (i=0; i<ndims; i++) {
        if (i != n) {
          ip[tick]    = mpi->ip[i];
          iproc[tick] = mpi->iproc[i];
          tick++;
        }
      }
      color = ArrayIndex1D(ndims-1,iproc,ip,NULL,0); 
      key   = mpi->ip[n];
      MPI_Comm_split(mpi->world,color,key,&mpi->comm[n]);
    }
    free(ip);
    free(iproc);
  }
  return(0);
}

int MPIFreeCommunicators(int ndims,void *m)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          n;
  
  for (n=0; n<ndims; n++) MPI_Comm_free(&mpi->comm[n]);
  free(mpi->comm);

  return(0);
}
