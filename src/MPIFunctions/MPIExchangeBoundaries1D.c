#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIExchangeBoundaries1D(void *m,double *x,int N,int ghosts,int dir,int ndims)
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) m;
  MPI_Request   *requests;
  MPI_Status    *statuses;
  int           ierr = 0, i;
  
  int *ip     = mpi->ip;
  int *iproc  = mpi->iproc;
  int non      = 0; /* number of neighbours */

  int neighbor_rank[2] = {-1,-1};
  int *nip             = (int*) calloc (ndims  ,sizeof(int));


  /* each process has 2 neighbors (except at physical boundaries)       */
  /* calculate the rank of these neighbors (-1 -> none)                 */
  ierr = ArrayCopy1D_int(ip,nip,ndims); CHECKERR(ierr); nip[dir]--;
  if (ip[dir] == 0)             neighbor_rank[0] = -1;
  else                          neighbor_rank[0] = MPIRank1D(ndims,iproc,nip);
  ierr = ArrayCopy1D_int(ip,nip,ndims); CHECKERR(ierr); nip[dir]++;
  if (ip[dir] == (iproc[dir]-1))neighbor_rank[1] = -1;
  else                          neighbor_rank[1] = MPIRank1D(ndims,iproc,nip);

  /* Allocate send and receive buffers */
  double **sendbuf,**recvbuf;
  sendbuf = (double**) calloc (2,sizeof(double*));
  recvbuf = (double**) calloc (2,sizeof(double*));
  sendbuf[0] = (double*) calloc(ghosts,sizeof(double));
  sendbuf[1] = (double*) calloc(ghosts,sizeof(double));
  recvbuf[0] = (double*) calloc(ghosts,sizeof(double));
  recvbuf[1] = (double*) calloc(ghosts,sizeof(double));

  /* count number of neighbors and copy data to send buffers */
  non = 0;
  if (neighbor_rank[0] != -1) {
    non++;
    for (i = 0; i < ghosts; i++) sendbuf[0][i] = x[i+ghosts];
  }
  if (neighbor_rank[1] != -1) {
    non++;
    for (i = 0; i < ghosts; i++) sendbuf[1][i] = x[i+N];
  }
  requests = (MPI_Request*) calloc(2*non,sizeof(MPI_Request));
  statuses = (MPI_Status* ) calloc(2*non,sizeof(MPI_Status ));

  /* exchange the data */
  int tick = 0;
  if (neighbor_rank[0]!= -1) {
    MPI_Irecv(recvbuf[0],ghosts,MPI_DOUBLE,neighbor_rank[0],1631,MPI_COMM_WORLD,&requests[tick]);
    MPI_Isend(sendbuf[0],ghosts,MPI_DOUBLE,neighbor_rank[0],1631,MPI_COMM_WORLD,&requests[tick+non]);
    tick++;
  }
  if (neighbor_rank[1] != -1) {
    MPI_Irecv(recvbuf[1],ghosts,MPI_DOUBLE,neighbor_rank[1],1631,MPI_COMM_WORLD,&requests[tick]);
    MPI_Isend(sendbuf[1],ghosts,MPI_DOUBLE,neighbor_rank[1],1631,MPI_COMM_WORLD,&requests[tick+non]);
    tick++;
  }

  /* Wait till data transfer is done */
  MPI_Waitall(2*non,requests,statuses);
  if (requests) free(requests);
  if (statuses) free(statuses);

  /* copy received data to ghost points */
  if (neighbor_rank[0] != -1) for (i = 0; i < ghosts; i++) x[i]          = recvbuf[0][i];
  if (neighbor_rank[1] != -1) for (i = 0; i < ghosts; i++) x[i+N+ghosts] = recvbuf[1][i];
  
  /* free send and receive buffers */
  free(sendbuf[0]);
  free(sendbuf[1]);
  free(recvbuf[0]);
  free(recvbuf[1]);
  free(sendbuf);
  free(recvbuf);

  /* free other allocated variables */
  free(nip);
#endif
  return(0);
}

