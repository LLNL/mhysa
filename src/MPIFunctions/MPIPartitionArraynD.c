#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIPartitionArraynD(int ndims,void *m,double *xg,double *x,int *dim_global,int *dim_local,
                        int ghosts,int nvars)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0;

  int *is, *ie, *index, *bounds;
  is     = (int*) calloc (ndims,sizeof(int));
  ie     = (int*) calloc (ndims,sizeof(int));
  index  = (int*) calloc (ndims,sizeof(int));
  bounds = (int*) calloc (ndims,sizeof(int));
  double *buffer;

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->rank);
    ierr = 1;
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array is not allocated on root processor.\n");
    ierr = 1;
  }

  if (!mpi->rank) {
    int proc;
    for (proc = 0; proc < mpi->nproc; proc++) {
      int d,done,size;
      /* Find out the domain limits for each process */
      ierr = MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie); CHECKERR(ierr);
      size = 1;
      for (d=0; d<ndims; d++) {
        size *= (ie[d]-is[d]);
        bounds[d] = ie[d] - is[d];
      }
      buffer = (double*) calloc (size*nvars,sizeof(double));
      done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
      while (!done) {
        int p1 = ArrayIndex1D(ndims,dim_global,index,is  ,0);
        int p2 = ArrayIndex1D(ndims,bounds    ,index,NULL,0);
        int v; for (v=0; v<nvars; v++) buffer[nvars*p2+v] = xg[nvars*p1+v];
        done = ArrayIncrementIndex(ndims,bounds,index);
      }
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size*nvars,MPI_DOUBLE,proc,1538,mpi->world);
#endif
      } else {
        done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
        while (!done) {
          int p1 = ArrayIndex1D(ndims,dim_local,index,NULL,ghosts);
          int p2 = ArrayIndex1D(ndims,dim_local,index,NULL,0     );
          int v; for (v=0; v<nvars; v++) x[nvars*p1+v] = buffer[nvars*p2+v];
          done = ArrayIncrementIndex(ndims,dim_local,index);
        }
      }
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Status  status;
    int d, done, size;
    size = 1; for (d=0; d<ndims; d++) size *= dim_local[d];
    buffer = (double*) calloc (size*nvars,sizeof(double));
    MPI_Recv(buffer,size*nvars,MPI_DOUBLE,0,1538,mpi->world,&status);
    done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
    while (!done) {
      int p1 = ArrayIndex1D(ndims,dim_local,index,NULL,ghosts);
      int p2 = ArrayIndex1D(ndims,dim_local,index,NULL,0     );
      int v; for (v=0; v<nvars; v++) x[nvars*p1+v] = buffer[nvars*p2+v];
      done = ArrayIncrementIndex(ndims,dim_local,index);
    }
    free(buffer);
#endif
  }
  free(is);
  free(ie);
  free(index);
  free(bounds);
  return(ierr);
}
