#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIGatherArray1D(void *m,double *xg,double *x,int istart,int iend,
                     int N_local,int ghosts)
{
  MPIVariables *mpi = (MPIVariables*) m;
  int          ierr = 0,i;
  double       *buffer;

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

  /* create and copy data to a buffer to send to the root process */
  buffer = (double*) calloc (N_local,sizeof(double));
  for (i=0; i<N_local; i++) buffer[i] = x[i+ghosts];

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
        for (i=0; i<size; i++) xg[is+i] = recvbuf[i];
        free(recvbuf);
#endif
      } else for (i=0; i<size; i++) xg[is+i] = buffer[i];
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
