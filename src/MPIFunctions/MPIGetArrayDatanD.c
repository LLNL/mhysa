#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIGetArrayDatanD(double *xbuf,double *x,int *source,int *dest,int *limits,
                      int *dim,int ghosts,int ndims,int nvars,void *m)
{
  MPIVariables *mpi  = (MPIVariables*) m;
  int          d;

  int source_rank = MPIRank1D(ndims,mpi->iproc,source);
  int dest_rank   = MPIRank1D(ndims,mpi->iproc,dest  );

  int is[ndims], ie[ndims], index[ndims], bounds[ndims], size;
  size    = 1;
  for (d=0; d<ndims; d++) {
    is[d] =  limits[2*d  ];
    ie[d] =  limits[2*d+1];
    size  *= (ie[d] - is[d]);
  }
  _ArraySubtract1D_(bounds,ie,is,ndims);

  if (source_rank == dest_rank) {
    /* source and dest are the same process */
    int done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p1; _ArrayIndex1D_(ndims,bounds,index,0,p1);
      int p2; _ArrayIndex1DWO_(ndims,dim,index,is,ghosts,p2);
      _ArrayCopy1D_((x+nvars*p2),(xbuf+nvars*p1),nvars);
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  } else {
#ifdef serial
    fprintf(stderr,"Error in MPIGetArrayDatanD(): This is a serial run. Source and ");
    fprintf(stderr,"destination ranks must be equal. Instead they are %d and %d.\n",
                    source_rank,dest_rank);
    return(1);
#else
    if (mpi->rank == source_rank) {
      double *buf = (double*) calloc (size*nvars, sizeof(double));
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1D_(ndims,bounds,index,0,p1);
        int p2; _ArrayIndex1DWO_(ndims,dim,index,is,ghosts,p2);
        _ArrayCopy1D_((x+nvars*p2),(buf+nvars*p1),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
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
  return(0);
}

