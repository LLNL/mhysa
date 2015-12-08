/*! @file MPIGatherArray1D.c
    @brief Gathers local 1D arrays to a global 1D array
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Gathers the contents of a 1D array (partitioned amongst MPI ranks) into a global 1D array
  on the root rank (rank 0). See documentation of MPIExchangeBoundaries1D() on what a "1D 
  array" is in the context of a multidimensional simulation. The 1D array must be the same 
  along spatial dimensions normal to the one it represents.

  Notes:
  + The global array must not have ghost points.
  + The global array must be preallocated on only rank 0. On other ranks, it must be NULL.
  + Since this function deals with a 1D array, more than one rank may be sending the same 
    piece of data to rank 0 (i.e. if there are more than one MPI rank along the dimensions
    normal to one corresponding to \a x ). The implementation of this function ignores this
    and overwrites that portion with the latest data sent.
*/
int MPIGatherArray1D(
                      void    *m,       /*!< MPI object of type #MPIVariables */
                      double  *xg,      /*!< Global 1D array (must be preallocated) without ghost points */
                      double  *x,       /*!< Local 1D array to be gathered */
                      int     istart,   /*!< Starting index (global) of this rank's portion of the array */
                      int     iend,     /*!< Ending index (global) of this rank's portion of the array + 1 */
                      int     N_local,  /*!< Local size of the array */
                      int     ghosts    /*!< Number of ghost points */
                    )
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0;

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIGatherArray1D(): global array exists on non-root processors (rank %d).\n",
            mpi->rank);
    ierr = 1;
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIGatherArray1D(): global array is not allocated on root processor.\n");
    ierr = 1;
  }

  /* create and copy data to a buffer to send to the root process */
  double *buffer = (double*) calloc (N_local,sizeof(double));
  _ArrayCopy1D_((x+ghosts),(buffer),N_local);

  if (!mpi->rank) {
#ifndef serial
    MPI_Status status;
#endif
    int proc;
    for (proc = 0; proc < mpi->nproc; proc++) {
      /* Find out the domain limits for each process */
      int is,ie;
      if (proc) {
#ifndef serial
        MPI_Recv(&is,1,MPI_INT,proc,1442,mpi->world,&status);
        MPI_Recv(&ie,1,MPI_INT,proc,1443,mpi->world,&status);
#endif
      } else { is = istart; ie = iend; }
      int size = ie - is;
      if (proc) {
#ifndef serial
        double *recvbuf = (double*) calloc (size,sizeof(double));
        MPI_Recv(recvbuf,size,MPI_DOUBLE,proc,1916,mpi->world,&status);
        _ArrayCopy1D_((recvbuf),(xg+is),size);
        free(recvbuf);
#endif
      } else _ArrayCopy1D_(buffer,(xg+is),size);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes - send stuff to root */
    MPI_Send(&istart,1      ,MPI_INT   ,0,1442,mpi->world);
    MPI_Send(&iend  ,1      ,MPI_INT   ,0,1443,mpi->world);
    MPI_Send(buffer ,N_local,MPI_DOUBLE,0,1916,mpi->world);
#endif
  }

  free(buffer);
  return(ierr);
}
