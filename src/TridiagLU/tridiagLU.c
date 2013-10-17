#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

int tridiagLU(double **a,double **b,double **c,double **x,
              int n,int ns,void *r,void *m)
{
  TridiagLUTime   *runtimes = (TridiagLUTime*) r;
  int             d,i,istart,iend;
  int             rank,nproc;
  struct timeval  start,stage1,stage2,stage3,stage4;
#ifndef serial
  MPIContext      *mpi = (MPIContext*) m;
  int             ierr = 0;
  const int       nvar = 4;
  int             *proc;
  double          *sendbuf,*recvbuf;
  MPI_Comm        *comm;
#endif

  /* start */
  gettimeofday(&start,NULL);

#ifdef serial
  rank  = 0;
  nproc = 1;
#else
  if (mpi) {
    rank  = mpi->rank;
    nproc = mpi->nproc;
    comm  = (MPI_Comm*) mpi->comm;
    proc  = mpi->proc;
    if (!proc) {
      fprintf(stderr,"Error in tridiagLU() on process %d: ",rank);
      fprintf(stderr,"aray \"proc\" is NULL.\n");
      return(-1);
    }
  } else {
    rank  = 0;
    nproc = 1;
    comm  = NULL;
    proc  = NULL;
  }
#endif


  if ((ns == 0) || (n == 0)) return(0);
  /* some allocations */
  double *xs1,*xp1; /* to exchange the first element on the next process   */
  xp1 = (double*) calloc(ns,sizeof(double)); for (d=0; d<ns; d++) xp1[d] = 0.0;
  xs1 = (double*) calloc(ns,sizeof(double)); for (d=0; d<ns; d++) xs1[d] = 0.0;

  /* Stage 1 - Parallel elimination of subdiagonal entries */
  istart  = (rank == 0 ? 1 : 2);
  iend    = n;
  for (d = 0; d < ns; d++) {
    for (i = istart; i < iend; i++) {
      if (b[d][i-1] == 0) return(-1);
      double factor = a[d][i] / b[d][i-1];
      b[d][i] -=  factor * c[d][i-1];
      a[d][i]  = -factor * a[d][i-1];
      x[d][i] -=  factor * x[d][i-1];
      if (rank) {
        double factor = c[d][0] / b[d][i-1];
        c[d][0]  = -factor * c[d][i-1];
        b[d][0] -=  factor * a[d][i-1];
        x[d][0] -=  factor * x[d][i-1];
      }
    }
  }

  /* end of stage 1 */
  gettimeofday(&stage1,NULL);

  /* Stage 2 - Eliminate the first sub- & super-diagonal entries */
  /* This needs the last (a,b,c,x) from the previous process     */
#ifndef serial
  sendbuf = (double*) calloc (ns*nvar,sizeof(double*));
  recvbuf = (double*) calloc (ns*nvar,sizeof(double*));
  for (d=0; d<ns; d++) {
    sendbuf[d*nvar+0] = a[d][n-1]; 
    sendbuf[d*nvar+1] = b[d][n-1]; 
    sendbuf[d*nvar+2] = c[d][n-1]; 
    sendbuf[d*nvar+3] = x[d][n-1];
  }
  if (nproc > 1) {
    MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if (rank)             MPI_Irecv(recvbuf,nvar*ns,MPI_DOUBLE,proc[rank-1],1436,*comm,&req[0]);
    if (rank != nproc-1)  MPI_Isend(sendbuf,nvar*ns,MPI_DOUBLE,proc[rank+1],1436,*comm,&req[1]);
    MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  }
  /* The first process sits this one out */
  if (rank) {
    for (d = 0; d < ns; d++) {
      double am1, bm1, cm1, xm1;
      am1 = recvbuf[d*nvar+0]; 
      bm1 = recvbuf[d*nvar+1]; 
      cm1 = recvbuf[d*nvar+2]; 
      xm1 = recvbuf[d*nvar+3];
      double factor;
      if (bm1 == 0) return(-1);
      factor =  a[d][0] / bm1;
      b[d][0]  -=  factor * cm1;
      a[d][0]   = -factor * am1;
      x[d][0]  -=  factor * xm1;
      if (b[d][n-1] == 0) return(-1);
      factor =  c[d][0] / b[d][n-1];
      b[d][0]  -=  factor * a[d][n-1];
      c[d][0]   = -factor * c[d][n-1];
      x[d][0]  -=  factor * x[d][n-1];
    }
  }
  free(sendbuf); free(recvbuf);

  /* end of stage 2 */
  gettimeofday(&stage2,NULL);

  /* Stage 3 - Solve the reduced (nproc-1) X (nproc-1) tridiagonal system   */
  if (nproc > 1) {
    double **zero, **one;
    zero    = (double**) calloc (ns,sizeof(double*));
    one     = (double**) calloc (ns,sizeof(double*));
    for (d=0; d<ns; d++) {
      zero[d] = (double* ) calloc (1,sizeof(double )); zero[d][0] = 0.0;
      one [d] = (double* ) calloc (1,sizeof(double )); one [d][0] = 1.0;
    }
    /* Solving the reduced system by gather-and-solve algorithm */
    if (rank) ierr = tridiagLUGS(a,b,c,x,1,ns,NULL,mpi);
    else      ierr = tridiagLUGS(zero,one,zero,zero,1,ns,NULL,mpi);
    if (ierr) return(ierr);
    for (d=0; d<ns; d++) free(zero[d]); free(zero);
    for (d=0; d<ns; d++) free(one [d]); free(one );

    /* Each process, get the first x of the next process */
    MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    for (d=0; d<ns; d++)  xs1[d] = x[d][0];
    if (rank+1 < nproc) MPI_Irecv(xp1,ns,MPI_DOUBLE,proc[rank+1],1323,*comm,&req[0]);
    if (rank)           MPI_Isend(xs1,ns,MPI_DOUBLE,proc[rank-1],1323,*comm,&req[1]);
    MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  }
#else
  if (nproc > 1) {
    fprintf(stderr,"Error: nproc > 1 for a serial run!\n");
    return(1);
  }
#endif /* if not serial */
  /* end of stage 3 */
  gettimeofday(&stage3,NULL);

  /* Stage 4 - Parallel back-substitution to get the solution  */
  istart = n-1;
  iend   = (rank == 0 ? 0 : 1);
  for (d = 0; d < ns; d++) {
    if (b[d][istart] == 0) return(-1);
    x[d][istart] = (x[d][istart]-a[d][istart]*x[d][0]-c[d][istart]*xp1[d]) / b[d][istart];
    for (i = istart-1; i > iend-1; i--) {
      if (b[d][i] == 0) return(-1);
      x[d][i] = (x[d][i]-c[d][i]*x[d][i+1]-a[d][i]*x[d][0]) / b[d][i];
    }
  }

  /* end of stage 4 */
  gettimeofday(&stage4,NULL);

  /* Done - now x contains the solution */
  free(xp1);
  free(xs1);

  /* save runtimes if needed */
  if (runtimes) {
    long long walltime;
    walltime = ((stage1.tv_sec * 1000000 + stage1.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
    runtimes->stage1_time = (double) walltime / 1000000.0;
    walltime = ((stage2.tv_sec * 1000000 + stage2.tv_usec) - (stage1.tv_sec * 1000000 + stage1.tv_usec));
    runtimes->stage2_time = (double) walltime / 1000000.0;
    walltime = ((stage3.tv_sec * 1000000 + stage3.tv_usec) - (stage2.tv_sec * 1000000 + stage2.tv_usec));
    runtimes->stage3_time = (double) walltime / 1000000.0;
    walltime = ((stage4.tv_sec * 1000000 + stage4.tv_usec) - (stage3.tv_sec * 1000000 + stage3.tv_usec));
    runtimes->stage4_time = (double) walltime / 1000000.0;
    walltime = ((stage4.tv_sec * 1000000 + stage4.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
    runtimes->total_time = (double) walltime / 1000000.0;
  }
  return(0);
}
