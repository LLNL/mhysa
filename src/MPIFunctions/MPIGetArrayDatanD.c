/*! @file MPIGetArrayDatanD.c
    @brief Send a part of a local array to another MPI rank
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  This function lets one rank get a portion of a local n-dimensional array on another rank. The n-dimensional
  array must be stored in the memory as a single-index array as described in the documentation of MPIExchangeBoundariesnD().
  The \a source rank sends to the \a dest rank a logically rectangular n-dimensional portion of its local copy of
  an array \a x. The extent of this logically rectangular portion is defined by \a limits.
  + \a limits is an array of size 2x the number of spatial dimensions, with elements as: 
    [ is[0], ie[0], is[1], ie[1], ..., is[ndims-1], ie[ndims-1] ], where is[n] is the starting index along 
    spatial dimension n, and ie[n] is the end index (+1) along spatial dimension n.
*/
int MPIGetArrayDatanD(
                        double  *xbuf,    /*!< preallocated memory on destination rank to hold the received data */
                        double  *x,       /*!< local array of which a part is needed */
                        int     *source,  /*!< MPI rank of the source */
                        int     *dest,    /*!< MPI rank of the destination */
                        int     *limits,  /*!< Integer array (of size 2*ndims) with the start and end indices 
                                               along each spatial dimension of the desired portion of the array */
                        int     *dim,     /*!< Integer array whose elements are the local size of x in each spatial dimension */
                        int     ghosts,   /*!< Number of ghost points */
                        int     ndims,    /*!< Number of spatial dimensions */
                        int     nvars,    /*!< Number of variables (vector components) */
                        void    *m        /*!< MPI object of type #MPIVariables */
                     )
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

