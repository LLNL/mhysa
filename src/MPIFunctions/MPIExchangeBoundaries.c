#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIExchangeBoundaries(int ndims,int nvars,int *dim,int ghosts,void *m,double *var)
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) m;
  MPI_Request   *requests;
  MPI_Status    *statuses;
  int           ierr = 0, d;
  
  int *ip     = mpi->ip;
  int *iproc  = mpi->iproc;
  int n_neighbors = 0;

  int *neighbor_rank  = (int*) calloc (2*ndims,sizeof(int));
  int *nip            = (int*) calloc (ndims  ,sizeof(int));
  int *index          = (int*) calloc (ndims  ,sizeof(int));
  int *bounds         = (int*) calloc (ndims,sizeof(int));
  int *offset         = (int*) calloc (ndims,sizeof(int));

  /* each process has 2*ndims neighbors (except at physical boundaries) */
  /* calculate the rank of these neighbors (-1 -> none)                 */
  for (d = 0; d < ndims; d++) {
    ierr = ArrayCopy1D_int(ip,nip,ndims); CHECKERR(ierr); nip[d]--;
    if (ip[d] == 0)             neighbor_rank[2*d]   = -1;
    else                        neighbor_rank[2*d]   = MPIRank1D(ndims,iproc,nip);
    ierr = ArrayCopy1D_int(ip,nip,ndims); CHECKERR(ierr); nip[d]++;
    if (ip[d] == (iproc[d]-1))  neighbor_rank[2*d+1] = -1;
    else                        neighbor_rank[2*d+1] = MPIRank1D(ndims,iproc,nip);
  }

  /* calculate dimensions of each of the send-receive regions */
  int *bufdim = (int*) calloc (ndims,sizeof(int));
  for (d = 0; d < ndims; d++) {
    bufdim[d] = 1;
    int i;
    for (i = 0; i < ndims; i++) {
      if (i == d) bufdim[d] *= ghosts;
      else        bufdim[d] *= dim[i];
    }
  }
  
  /* Allocate send and receive buffers */
  double **sendbuf,**recvbuf;
  sendbuf = (double**) calloc (2*ndims,sizeof(double*));
  recvbuf = (double**) calloc (2*ndims,sizeof(double*));
  for (d = 0; d < ndims; d++) {
    sendbuf[2*d]   = (double*) calloc(bufdim[d]*nvars,sizeof(double));
    sendbuf[2*d+1] = (double*) calloc(bufdim[d]*nvars,sizeof(double));
    recvbuf[2*d]   = (double*) calloc(bufdim[d]*nvars,sizeof(double));
    recvbuf[2*d+1] = (double*) calloc(bufdim[d]*nvars,sizeof(double));
  }

  /* count number of neighbors and copy data to send buffers */
  n_neighbors = 0;
  for (d = 0; d < ndims; d++) {
    ierr = ArrayCopy1D_int(dim,bounds,ndims); CHECKERR(ierr); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      n_neighbors++;
      ierr = ArraySetValue_int(offset,ndims,0); CHECKERR(ierr);
      int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
      while (!done) {
        int p1 = ArrayIndex1D(ndims,dim   ,index,offset,ghosts);
        int p2 = ArrayIndex1D(ndims,bounds,index,NULL  ,0     );
        int v; for (v=0; v<nvars; v++) sendbuf[2*d][nvars*p2+v] = var[nvars*p1+v];
        done = ArrayIncrementIndex(ndims,bounds,index);
      }
    }
    if (neighbor_rank[2*d+1] != -1) {
      n_neighbors++;
      ierr = ArraySetValue_int(offset,ndims,0); CHECKERR(ierr); offset[d] = dim[d]-ghosts;
      int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
      while (!done) {
        int p1 = ArrayIndex1D(ndims,dim   ,index,offset,ghosts);
        int p2 = ArrayIndex1D(ndims,bounds,index,NULL  ,0     );
        int v; for (v=0; v<nvars; v++) sendbuf[2*d+1][nvars*p2+v] = var[nvars*p1+v];
        done = ArrayIncrementIndex(ndims,bounds,index);
      }
    }
  }
  requests = (MPI_Request*) calloc(2*n_neighbors,sizeof(MPI_Request));
  statuses = (MPI_Status* ) calloc(2*n_neighbors,sizeof(MPI_Status ));

  /* exchange the data */
  int tick = 0;
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Irecv(recvbuf[2*d  ],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                MPI_COMM_WORLD,&requests[tick]);
      MPI_Isend(sendbuf[2*d  ],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                MPI_COMM_WORLD,&requests[tick+n_neighbors]);
      tick++;
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Irecv(recvbuf[2*d+1],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                MPI_COMM_WORLD,&requests[tick]);
      MPI_Isend(sendbuf[2*d+1],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                MPI_COMM_WORLD,&requests[tick+n_neighbors]);
      tick++;
    }
  }

  /* Wait till data transfer is done */
  MPI_Waitall(2*n_neighbors,requests,statuses);
  if (requests) free(requests);
  if (statuses) free(statuses);

  /* copy received data to ghost points */
  for (d = 0; d < ndims; d++) {
    ierr = ArrayCopy1D_int(dim,bounds,ndims); CHECKERR(ierr); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      ierr = ArraySetValue_int(offset,ndims,0); CHECKERR(ierr); offset[d] = -ghosts;
      int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
      while (!done) {
        int p1 = ArrayIndex1D(ndims,dim   ,index,offset,ghosts);
        int p2 = ArrayIndex1D(ndims,bounds,index,NULL  ,0     );
        int v; for (v=0; v<nvars; v++) var[nvars*p1+v] = recvbuf[2*d][nvars*p2+v];
        done = ArrayIncrementIndex(ndims,bounds,index);
      }
    }
    if (neighbor_rank[2*d+1] != -1) {
      ierr = ArraySetValue_int(offset,ndims,0); CHECKERR(ierr); offset[d] = dim[d];
      int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
      while (!done) {
        int p1 = ArrayIndex1D(ndims,dim   ,index,offset,ghosts);
        int p2 = ArrayIndex1D(ndims,bounds,index,NULL  ,0     );
        int v; for (v=0; v<nvars; v++) var[nvars*p1+v] = recvbuf[2*d+1][nvars*p2+v];
        done = ArrayIncrementIndex(ndims,bounds,index);
      }
    }
  }
  
  /* free send and receive buffers */
  for (d = 0; d < ndims; d++) {
    free(sendbuf[2*d  ]);
    free(sendbuf[2*d+1]);
    free(recvbuf[2*d  ]);
    free(recvbuf[2*d+1]);
  }
  free(sendbuf);
  free(recvbuf);

  /* free other allocated variables */
  free(bounds);
  free(offset);
  free(index);
  free(nip);
  free(neighbor_rank);
#endif
  return(0);
}

