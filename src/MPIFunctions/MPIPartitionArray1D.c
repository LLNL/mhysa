#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIPartitionArray1D(void *m,int nproc,int rank,double *xg,double *x,
                        int N_global,int N_local,int ghosts_global,int ghosts_local)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0;
  double       *buffer;

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array exists on non-root processors.\n");
    ierr = 1;
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array is not allocated on root processor.\n");
    ierr = 1;
  }

  if (!mpi->rank) {
    int     is, ie, proc;
    is = ie = 0;
    for (proc = 0; proc < mpi->nproc; proc++) {
      int i,size;
      /* Find out the domain limits for each process */
      size = MPIPartition1D(N_global,nproc,proc);
      ie = is + size;
      buffer = (double*) calloc (size,sizeof(double));
      for (i=0; i<size; i++) buffer[i] = xg[i+is+ghosts_global];
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size,MPI_DOUBLE,proc,1538,MPI_COMM_WORLD);
#else
        fprintf(stderr,"Error in MPIPartitionArray3D(): This is a serial run; number of processes ");
        fprintf(stderr,"should be 1. Instead it is %d.\n",mpi->nproc);
        return(1);
#endif
      } else for (i=0; i<N_local; i++) x[i] = buffer[i];
      free(buffer);
      is += size;
    }

  } else {
    int i;
#ifdef serial
    fprintf(stderr,"Error in MPIPartitionArray3D(): This is a serial run; process rank ");
    fprintf(stderr,"should be 0. Instead it is %d.\n",mpi->rank);
    return(1);
#else
    /* Meanwhile, on other processes */
    MPI_Status  status;
    buffer = (double*) calloc (N_local,sizeof(double));
    MPI_Recv(buffer,N_local,MPI_DOUBLE,0,1538,MPI_COMM_WORLD,&status);
    for (i=0; i<N_local; i++) x[i+ghosts_local] = buffer[i];
    free(buffer);
#endif
  }
  return(ierr);
}
