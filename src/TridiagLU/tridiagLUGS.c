#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

int tridiagLUGS(double **a,double **b,double **c,double **x,
              int n,int ns,void *r,void *m)
{
#ifdef serial

  /* Serial compilation */
  return(tridiagLU(a,b,c,x,n,ns,r,m));

#else

  /* Parallel compilation */
  MPIContext      *mpi = (MPIContext*) m;

  /* MPIContext argument = NULL => Serial solve */
  if (!mpi) return(tridiagLU(a,b,c,x,n,ns,r,NULL));

  int         d,i,ierr = 0,dstart,istart,p,q;
  const int   nvar = 4;
  double      *sendbuf,*recvbuf;
  int         rank  = mpi->rank;
  int         nproc = mpi->nproc;
  int         *proc = mpi->proc;
  MPI_Comm    *comm = (MPI_Comm*) mpi->comm;

  if (!proc) {
    fprintf(stderr,"Error in tridiagLUGS() on process %d: ",rank);
    fprintf(stderr,"array \"proc\" is NULL.\n");
    return(-1);
  }

  if ((ns == 0) || (n == 0)) return(0);

  /* 
    each process needs to know the local sizes of every other process 
    and total size of the system
  */
  int *N = (int*) calloc (nproc, sizeof(int));
  MPI_Allgather(&n,1,MPI_INT,N,1,MPI_INT,*comm);
  int NT = 0; for (i=0; i<nproc; i++) NT += N[i];

  /* counts and displacements for gather and scatter operations */
  int *counts = (int*) calloc (nproc, sizeof(int));
  int *displ  = (int*) calloc (nproc, sizeof(int));

  /* on all processes, calculate the number of systems each     */
  /* process has to solve                                       */
  int *ns_local = (int*) calloc (nproc,sizeof(int));
  for (p=0; p<nproc; p++)    ns_local[p] = ns / nproc; 
  for (p=0; p<ns%nproc; p++) ns_local[p]++;

  /* allocate the arrays for the gathered tridiagonal systems */
  double **ra,**rb,**rc,**rx; 
  if (ns_local[rank] > 0) {
    ra = (double**) calloc (ns_local[rank],sizeof(double));
    rb = (double**) calloc (ns_local[rank],sizeof(double));
    rc = (double**) calloc (ns_local[rank],sizeof(double));
    rx = (double**) calloc (ns_local[rank],sizeof(double));
    for (d = 0; d < ns_local[rank]; d++) {
      ra[d] = (double*) calloc (NT, sizeof(double));
      rb[d] = (double*) calloc (NT, sizeof(double));
      rc[d] = (double*) calloc (NT, sizeof(double));
      rx[d] = (double*) calloc (NT, sizeof(double));
      for (i = 0; i < NT; i++) 
        ra[d][i] = rb[d][i] = rc[d][i] = rx[d][i] = 0.0;
    }
  }

  /* Gather the systems on each process */
  /* allocate receive buffer */
  if (ns_local[rank] > 0) 
    recvbuf = (double*) calloc (ns_local[rank]*nvar*NT,sizeof(double));
  else recvbuf = NULL;
  dstart = 0;
  for (p = 0; p < nproc; p++) {
    if (ns_local[p] > 0) {
      /* allocate send buffer and form the send packet of data */
      sendbuf = (double*) calloc (nvar*n*ns_local[p],sizeof(double));
      for (d = 0; d < ns_local[p]; d++) {
        for (i = 0; i < n; i++) {
          sendbuf[n*nvar*d+n*0+i] = a[d+dstart][i];
          sendbuf[n*nvar*d+n*1+i] = b[d+dstart][i];
          sendbuf[n*nvar*d+n*2+i] = c[d+dstart][i];
          sendbuf[n*nvar*d+n*3+i] = x[d+dstart][i];
        }
      }
      dstart += ns_local[p];

      /* gather these reduced systems on process with rank = p */
      if (rank == p) {
        for (q = 0; q < nproc; q++) {
          counts[q] = nvar*N[q]*ns_local[p];
          displ [q] = (q == 0 ? 0 : displ[q-1]+counts[q-1]);
        }
      }
      MPI_Gatherv(sendbuf,nvar*n*ns_local[p],MPI_DOUBLE,
                  recvbuf,counts,displ,MPI_DOUBLE,proc[p],*comm);

      /* deallocate send buffer */
      free(sendbuf);
    }
  }
  /* extract the data from the recvbuf and solve */
  istart = 0;
  for (q = 0; q < nproc; q++) {
    for (d = 0; d < ns_local[rank]; d++) {
      for (i = 0; i < N[q]; i++) {
        ra[d][istart+i] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+0*N[q]+i];
        rb[d][istart+i] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+1*N[q]+i];
        rc[d][istart+i] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+2*N[q]+i];
        rx[d][istart+i] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+3*N[q]+i];
      }
    }
    istart += N[q];
  }
  /* deallocate receive buffer */
  if (recvbuf)  free(recvbuf);

  /* solve the gathered systems in serial */
  ierr = tridiagLU(ra,rb,rc,rx,NT,ns_local[rank],NULL,NULL);
  if (ierr) return(ierr);

  /* allocate send buffer and save the data to send */
  if (ns_local[rank] > 0)
    sendbuf = (double*) calloc (ns_local[rank]*NT,sizeof(double));
  else sendbuf = NULL;
  istart = 0;
  for (q = 0; q < nproc; q++) {
    for (i = 0; i < N[q]; i++) {
      for (d = 0; d < ns_local[rank]; d++) {
        sendbuf[istart*ns_local[rank]+d*N[q]+i] = rx[d][istart+i];
      }
    }
    istart += N[q];
  }
  dstart = 0;
  for (p = 0; p < nproc; p++) {
    if (ns_local[p] > 0) {

      /* allocate receive buffer */
      recvbuf = (double*) calloc (ns_local[p]*n, sizeof(double));

      /* scatter the solution back */
      for (q = 0; q < nproc; q++) {
        counts[q] = ns_local[p]*N[q];
        displ[q]  = (q == 0 ? 0 : displ[q-1]+counts[q-1]);
      }
      MPI_Scatterv(sendbuf,counts,displ,MPI_DOUBLE,
                   recvbuf,ns_local[p]*n,MPI_DOUBLE,
                   proc[p],*comm);
      /* save the solution on all root processes */
      for (d = 0; d < ns_local[p]; d++) {
        for (i = 0; i < n; i++) {
          x[d+dstart][i] = recvbuf[d*n+i];
        }
      }
      dstart += ns_local[p];
      /* deallocate receive buffer */
      free(recvbuf);
    }
  }
  /* deallocate send buffer */
  if (sendbuf) free(sendbuf);

  /* clean up */
  if (ns_local[rank] > 0) {
    for (d = 0; d < ns_local[rank]; d++) {
      free(ra[d]);
      free(rb[d]);
      free(rc[d]);
      free(rx[d]);
    }
    free(ra);
    free(rb);
    free(rc);
    free(rx);
  }
  free(ns_local);
  free(N);
  free(displ);
  free(counts);

  return(0);
#endif
}
