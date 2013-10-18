#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIGetArrayDatanD(double *xbuf,double *x,int *source,int *dest,int *limits,
                      int *dim,int ghosts,int ndims,int nvars,void *m)
{
  MPIVariables *mpi  = (MPIVariables*) m;
  int          ierr = 0, d;

  int source_rank = MPIRank1D(ndims,mpi->iproc,source);
  int dest_rank   = MPIRank1D(ndims,mpi->iproc,dest  );

  int *is, *ie, *bounds, *index, size;
  is      = (int*) calloc (ndims,sizeof(int));
  ie      = (int*) calloc (ndims,sizeof(int));
  bounds  = (int*) calloc (ndims,sizeof(int));
  index   = (int*) calloc (ndims,sizeof(int));
  size    = 1;
  for (d=0; d<ndims; d++) {
    is[d] =  limits[2*d  ];
    ie[d] =  limits[2*d+1];
    size  *= (ie[d] - is[d]);
  }
  ierr = ArraySubtract1D_int(bounds,ie,is,ndims);

  if (source_rank == dest_rank) {
    /* source and dest are the same process */
    int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
    while (!done) {
      int p1 = ArrayIndex1D(ndims,bounds,index,NULL,0     );
      int p2 = ArrayIndex1D(ndims,dim   ,index,is  ,ghosts);
      int v; for (v=0; v<nvars; v++) xbuf[nvars*p1+v] = x[nvars*p2+v];
      done = ArrayIncrementIndex(ndims,bounds,index);
    }
  } else {
#ifdef serial
    fprintf(stderr,"Error in MPIGetArrayDatanD(): This is a serial run. Source and ");
    fprintf(stderr,"destination ranks must be equal. Instead they are %d and %d.\n",
                    source_rank,dest_rank);
    return(1);
#else
    if (mpi->rank == source_rank) {
      double *buf = (double*) calloc (size*nvars,sizeof(double));
      int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
      while (!done) {
        int p1 = ArrayIndex1D(ndims,bounds,index,NULL,0     );
        int p2 = ArrayIndex1D(ndims,dim   ,index,is  ,ghosts);
        int v; for (v=0; v<nvars; v++) buf[nvars*p1+v] = x[nvars*p2+v];
        done = ArrayIncrementIndex(ndims,bounds,index);
      }
      MPI_Send(buf,size*nvars,MPI_DOUBLE,dest_rank,2211,mpi->world);
      free(buf);
    } else if (mpi->rank == dest_rank) {
      MPI_Status status;
      MPI_Recv(xbuf,size*nvars,MPI_DOUBLE,source_rank,2211,mpi->world,&status);
    } else {
      fprintf(stderr,"Error in MPIGetArrayData3D(): Process %d shouldn't have entered this function.\n",
              mpi->rank);
      return(1);
    }
#endif
  }

  free(is);
  free(ie);
  free(bounds);
  free(index);
  return(0);
}

