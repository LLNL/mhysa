/*! @file MPIGatherArraynD.c
    @brief Gather an n-dimensional array in to a global array
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Gathers the contents of an n-dimensional array, partitioned amongst the MPI ranks, in to a global
  array on rank 0. See documentation of MPIExchangeBoundariesnD() for how the n-dimensional array is
  stored in the memory as a single-index array.

  Notes:
  + The global array must have no ghost points.
  + The global array must be allocated only on rank 0. On other ranks, it must be NULL.
*/
int MPIGatherArraynD(
                      int     ndims,        /*!< Number of spatial dimensions */
                      void    *m,           /*!< MPI object of type #MPIVariables */
                      double  *xg,          /*!< Global array (preallocated) without ghost points */
                      double  *x,           /*!< Local array */
                      int     *dim_global,  /*!< Integer array with elements as global size along each spatial dimension */
                      int     *dim_local,   /*!< Integer array with elements as local size along each spatial dimension */
                      int     ghosts,       /*!< Number of ghost points */
                      int     nvars         /*!< Number of variables (vector components) */
                    )
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          d, size;
  _DECLARE_IERR_;

  int is[ndims], ie[ndims], index[ndims], bounds[ndims];

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIGatherArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->rank);
    return(1);
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIGatherArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

  /* calculate total size of local domain (w/o ghosts) */
  size = 1;
  for (d = 0; d < ndims; d++) size *= dim_local[d];

  /* create and copy data to send to root process */
  double *buffer = (double*) calloc (size*nvars, sizeof(double));
  IERR ArrayCopynD(ndims,x,buffer,dim_local,ghosts,0,index,nvars); CHECKERR(ierr);

  if (!mpi->rank) {
    int proc;
    for (proc = 0; proc < mpi->nproc; proc++) {
      int d,done,size;
      /* Find out the domain limits for each process */
      IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie); CHECKERR(ierr);
      size = 1;
      for (d=0; d<ndims; d++) {
        size *= (ie[d]-is[d]);
        bounds[d] = ie[d] - is[d];
      }
      if (proc) {
#ifndef serial
        MPI_Status status;
        double *recvbuf = (double*) calloc (size*nvars, sizeof(double));
        MPI_Recv(recvbuf,size*nvars,MPI_DOUBLE,proc,1902,mpi->world,&status);
        int done = 0; _ArraySetValue_(index,ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(ndims,bounds,index,0,p1);
          int p2; _ArrayIndex1DWO_(ndims,dim_global,index,is,0,p2);
          _ArrayCopy1D_((recvbuf+nvars*p1),(xg+nvars*p2),nvars);
          _ArrayIncrementIndex_(ndims,bounds,index,done);
        }
        free(recvbuf);
#endif
      } else {
        done = 0; _ArraySetValue_(index,ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(ndims,bounds,index,0,p1);
          int p2; _ArrayIndex1DWO_(ndims,dim_global,index,is,0,p2);
          _ArrayCopy1D_((buffer+nvars*p1),(xg+nvars*p2),nvars);
          _ArrayIncrementIndex_(ndims,bounds,index,done);
        }
      }
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Send(buffer,size*nvars,MPI_DOUBLE,0,1902,mpi->world);
#endif
  }

  free(buffer);
  return(0);
}
