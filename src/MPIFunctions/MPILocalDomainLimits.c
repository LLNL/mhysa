#include <stdlib.h>
#include <mpivars.h>

int MPILocalDomainLimits(int ndims,int p,void *m,int *dim_global,int *is, int *ie) 
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0,i;

  int *ip = (int*) calloc (ndims,sizeof(int));
  ierr = MPIRanknD(ndims,p,mpi->iproc,ip); if (ierr) return(ierr);

  for (i=0; i<ndims; i++) {
    int imax_local, isize, root = 0;
    ierr = MPIPartition1D(1,&dim_global[i],&mpi->iproc[i],&root ,&imax_local); if (ierr) return(ierr);
    ierr = MPIPartition1D(1,&dim_global[i],&mpi->iproc[i],&ip[i],&isize);      if (ierr) return(ierr);
    if (is)  is[i] = ip[i]*imax_local;
    if (ie)  ie[i] = ip[i]*imax_local + isize;
  }
  free(ip);
  return(0);
}
