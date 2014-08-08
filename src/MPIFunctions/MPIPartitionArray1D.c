#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIPartitionArray1D(void *m,double *xg,double *x,int istart,int iend,
                        int N_local,int ghosts)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0,i;
#ifndef serial
  MPI_Status   status;
#endif

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIPartitionArray1D(): global array exists on non-root processors (rank %d).\n",
            mpi->rank);
    ierr = 1;
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIPartitionArray1D(): global array is not allocated on root processor.\n");
    ierr = 1;
  }

  if (!mpi->rank) {
    int proc;
    for (proc = 0; proc < mpi->nproc; proc++) {
      /* Find out the domain limits for each process */
      int is,ie;
      if (proc) {
#ifndef serial
        MPI_Recv(&is,1,MPI_INT,proc,1442,mpi->world,&status);
        MPI_Recv(&ie,1,MPI_INT,proc,1443,mpi->world,&status);
#endif
      } else {
        is = istart;
        ie = iend;
      }
      int size = ie - is;
      double *buffer = (double*) calloc (size, sizeof(double));
      _ArrayCopy1D_((xg+is),buffer,size);
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size,MPI_DOUBLE,proc,1539,mpi->world);
#endif
      } else _ArrayCopy1D_(buffer,x,N_local);
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    /* send local start and end indices to root */
    MPI_Send(&istart,1,MPI_INT,0,1442,mpi->world);
    MPI_Send(&iend  ,1,MPI_INT,0,1443,mpi->world);
    double *buffer = (double*) calloc (N_local, sizeof(buffer));
    MPI_Recv(buffer,N_local,MPI_DOUBLE,0,1539,mpi->world,&status);
    _ArrayCopy1D_(buffer,x,N_local);
    free(buffer);
#endif
  }
  return(ierr);
}
