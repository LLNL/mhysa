#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIPartitionArray1D(void *m,double *xg,double *x,int istart,int iend,
                        int N_local,int ghosts_global,int ghosts_local)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0,i;
  double       *buffer;
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
        MPI_Recv(&is,1,MPI_INT,proc,1442,MPI_COMM_WORLD,&status);
        MPI_Recv(&ie,1,MPI_INT,proc,1443,MPI_COMM_WORLD,&status);
      } else {
        is = istart;
        ie = iend;
      }
      int size = ie - is;
      buffer = (double*) calloc (size,sizeof(double));
      for (i=0; i<size; i++) buffer[i] = xg[i+is+ghosts_global];
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size,MPI_DOUBLE,proc,1539,MPI_COMM_WORLD);
#endif
      } else for (i=0; i<N_local; i++) x[i] = buffer[i];
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    /* send local start and end indices to root */
    MPI_Send(&istart,1,MPI_INT,0,1442,MPI_COMM_WORLD);
    MPI_Send(&iend  ,1,MPI_INT,0,1443,MPI_COMM_WORLD);
    buffer = (double*) calloc (N_local,sizeof(double));
    MPI_Recv(buffer,N_local,MPI_DOUBLE,0,1539,MPI_COMM_WORLD,&status);
    for (i=0; i<N_local; i++) x[i+ghosts_local] = buffer[i];
    free(buffer);
#endif
  }
  return(ierr);
}
