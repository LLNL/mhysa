#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>

int MPILocalDomainLimits(int ndims,int p,void *m,int *dim_global,int *is, int *ie) 
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0,i;

  int *ip = (int*) calloc (ndims,sizeof(int));
  ierr = MPIRanknD(ndims,p,mpi->iproc,ip); CHECKERR(ierr);

  for (i=0; i<ndims; i++) {
    int imax_local, isize, root = 0;
    imax_local = MPIPartition1D(dim_global[i],mpi->iproc[i],root );
    isize      = MPIPartition1D(dim_global[i],mpi->iproc[i],ip[i]);
    if (is)  is[i] = ip[i]*imax_local;
    if (ie)  ie[i] = ip[i]*imax_local + isize;
  }
  free(ip);
  return(0);
}
