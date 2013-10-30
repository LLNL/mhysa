#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

int MPIExchangeBoundariesnD(int ndims,int nvars,int *dim,int ghosts,void *m,double *var)
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) m;
  MPI_Request   *rcvreq,*sndreq;
  int           ierr = 0, d;
  
  int *ip     = mpi->ip;
  int *iproc  = mpi->iproc;
  int *bcflag = mpi->bcperiodic;

  int *neighbor_rank  = (int*) calloc (2*ndims,sizeof(int));
  int *nip            = (int*) calloc (ndims  ,sizeof(int));
  int *index          = (int*) calloc (ndims  ,sizeof(int));
  int *bounds         = (int*) calloc (ndims,sizeof(int));
  int *offset         = (int*) calloc (ndims,sizeof(int));

  rcvreq = (MPI_Request*) calloc(2*ndims,sizeof(MPI_Request));
  sndreq = (MPI_Request*) calloc(2*ndims,sizeof(MPI_Request));
  for (d=0; d<2*ndims; d++) rcvreq[d] = sndreq[d] = MPI_REQUEST_NULL;

  /* each process has 2*ndims neighbors (except at non-periodic physical boundaries)  */
  /* calculate the rank of these neighbors (-1 -> none)                               */
  for (d = 0; d < ndims; d++) {
    ierr = ArrayCopy1D_int(ip,nip,ndims); CHECKERR(ierr); 
    if (ip[d] == 0) nip[d] = iproc[d]-1;
    else            nip[d]--;
    if ((ip[d] == 0) && (!bcflag[d])) neighbor_rank[2*d]   = -1;
    else                              neighbor_rank[2*d]   = MPIRank1D(ndims,iproc,nip);
    ierr = ArrayCopy1D_int(ip,nip,ndims); CHECKERR(ierr); 
    if (ip[d] == (iproc[d]-1)) nip[d] = 0;
    else                       nip[d]++;
    if ((ip[d] == (iproc[d]-1)) && (!bcflag[d]))  neighbor_rank[2*d+1] = -1;
    else                                          neighbor_rank[2*d+1] = MPIRank1D(ndims,iproc,nip);
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

  /* post the receive requests */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Irecv(recvbuf[2*d  ],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                mpi->world,&rcvreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Irecv(recvbuf[2*d+1],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1631,
                mpi->world,&rcvreq[2*d+1]);
    }
  }

  /* count number of neighbors and copy data to send buffers */
  for (d = 0; d < ndims; d++) {
    ierr = ArrayCopy1D_int(dim,bounds,ndims); CHECKERR(ierr); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
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

  /* send the data */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Isend(sendbuf[2*d  ],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1631,
                mpi->world,&sndreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Isend(sendbuf[2*d+1],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                mpi->world,&sndreq[2*d+1]);
    }
  }

  /* Wait till data is done received */
  MPI_Waitall(2*ndims,rcvreq,MPI_STATUS_IGNORE);

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
  
  /* Wait till send requests are complete before freeing memory */
  MPI_Waitall(2*ndims,sndreq,MPI_STATUS_IGNORE);

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
  free(bufdim);
  free(bounds);
  free(offset);
  free(index);
  free(nip);
  free(neighbor_rank);
  free(rcvreq);
  free(sndreq);
#endif
  return(0);
}

