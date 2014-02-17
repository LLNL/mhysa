#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

int tridiagLU(double *a,double *b,double *c,double *x,
              int n,int ns,void *r,void *m)
{
  TridiagLU       *params = (TridiagLU*) r;
  int             d,i,istart,iend;
  int             rank,nproc;
  struct timeval  start,stage1,stage2,stage3,stage4;

#ifdef serial
  rank  = 0;
  nproc = 1;
#else
  MPI_Comm        *comm = (MPI_Comm*) m;
  int             ierr = 0;
  const int       nvar = 4;

  if (comm) {
    MPI_Comm_size(*comm,&nproc);
    MPI_Comm_rank(*comm,&rank);
  } else {
    rank  = 0;
    nproc = 1;
  }
#endif

  if (!params) {
    fprintf(stderr,"Error in tridiagLU(): NULL pointer passed for parameters.\n");
    return(1);
  }

  /* start */
  gettimeofday(&start,NULL);

  if ((ns == 0) || (n == 0)) return(0);
  double *xs1, *xp1;
  xs1 = (double*) calloc (ns, sizeof(double));
  xp1 = (double*) calloc (ns, sizeof(double));
  for (i=0; i<ns; i++) xs1[i] = xp1[i] = 0;

  /* Stage 1 - Parallel elimination of subdiagonal entries */
  istart  = (rank == 0 ? 1 : 2);
  iend    = n;
  for (i = istart; i < iend; i++) {
    for (d = 0; d < ns; d++) {
      if (b[(i-1)*ns+d] == 0) return(-1);
      double factor = a[i*ns+d] / b[(i-1)*ns+d];
      b[i*ns+d] -=  factor * c[(i-1)*ns+d];
      a[i*ns+d]  = -factor * a[(i-1)*ns+d];
      x[i*ns+d] -=  factor * x[(i-1)*ns+d];
      if (rank) {
        double factor = c[d] / b[(i-1)*ns+d];
        c[d]  = -factor * c[(i-1)*ns+d];
        b[d] -=  factor * a[(i-1)*ns+d];
        x[d] -=  factor * x[(i-1)*ns+d];
      }
    }
  }

  /* end of stage 1 */
  gettimeofday(&stage1,NULL);

  /* Stage 2 - Eliminate the first sub- & super-diagonal entries */
  /* This needs the last (a,b,c,x) from the previous process     */
#ifndef serial
  double *sendbuf, *recvbuf;
  sendbuf = (double*) calloc (ns*nvar, sizeof(double));
  recvbuf = (double*) calloc (ns*nvar, sizeof(double));
  for (d=0; d<ns; d++) {
    sendbuf[d*nvar+0] = a[(n-1)*ns+d]; 
    sendbuf[d*nvar+1] = b[(n-1)*ns+d]; 
    sendbuf[d*nvar+2] = c[(n-1)*ns+d]; 
    sendbuf[d*nvar+3] = x[(n-1)*ns+d];
  }
  if (nproc > 1) {
    MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if (rank)             MPI_Irecv(recvbuf,nvar*ns,MPI_DOUBLE,rank-1,1436,*comm,&req[0]);
    if (rank != nproc-1)  MPI_Isend(sendbuf,nvar*ns,MPI_DOUBLE,rank+1,1436,*comm,&req[1]);
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
      factor =  a[d] / bm1;
      b[d]  -=  factor * cm1;
      a[d]   = -factor * am1;
      x[d]  -=  factor * xm1;
      if (b[(n-1)*ns+d] == 0) return(-1);
      factor =  c[d] / b[(n-1)*ns+d];
      b[d]  -=  factor * a[(n-1)*ns+d];
      c[d]   = -factor * c[(n-1)*ns+d];
      x[d]  -=  factor * x[(n-1)*ns+d];
    }
  }
  free(sendbuf);
  free(recvbuf);
#endif

  /* end of stage 2 */
  gettimeofday(&stage2,NULL);

  /* Stage 3 - Solve the reduced (nproc-1) X (nproc-1) tridiagonal system   */
#ifndef serial
  if (nproc > 1) {
    double *zero, *one;
    zero = (double*) calloc (ns, sizeof(double));
    one  = (double*) calloc (ns, sizeof(double));
    for (d=0; d<ns; d++) {
      zero[d] = 0.0;
      one [d] = 1.0;
    }
    if (!strcmp(params->reducedsolvetype,_TRIDIAG_GS_)) {
      /* Solving the reduced system by gather-and-solve algorithm */
      if (rank) ierr = tridiagLUGS(a,b,c,x,1,ns,params,comm);
      else      ierr = tridiagLUGS(zero,one,zero,zero,1,ns,params,comm);
      if (ierr) return(ierr);
    } else if (!strcmp(params->reducedsolvetype,_TRIDIAG_JACOBI_)) {
      /* Solving the reduced system iteratively with the Jacobi method */
      if (rank) ierr = tridiagIterJacobi(a,b,c,x,1,ns,params,comm);
      else      ierr = tridiagIterJacobi(zero,one,zero,zero,1,ns,params,comm);
    }
    free(zero);
    free(one);

    /* Each process, get the first x of the next process */
    MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    for (d=0; d<ns; d++)  xs1[d] = x[d];
    if (rank+1 < nproc) MPI_Irecv(xp1,ns,MPI_DOUBLE,rank+1,1323,*comm,&req[0]);
    if (rank)           MPI_Isend(xs1,ns,MPI_DOUBLE,rank-1,1323,*comm,&req[1]);
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
    if (b[istart*ns+d] == 0) return(-1);
    x[istart*ns+d] = (x[istart*ns+d]-a[istart*ns+d]*x[d]-c[istart*ns+d]*xp1[d]) / b[istart*ns+d];
  }
  for (i = istart-1; i > iend-1; i--) {
    for (d = 0; d < ns; d++) {
      if (b[i*ns+d] == 0) return(-1);
      x[i*ns+d] = (x[i*ns+d]-c[i*ns+d]*x[(i+1)*ns+d]-a[i*ns+d]*x[d]) / b[i*ns+d];
    }
  }

  /* end of stage 4 */
  gettimeofday(&stage4,NULL);

  /* Done - now x contains the solution */
  free(xs1);
  free(xp1);

  /* save runtimes if needed */
  long long walltime;
  walltime = ((stage1.tv_sec * 1000000 + stage1.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
  params->stage1_time = (double) walltime / 1000000.0;
  walltime = ((stage2.tv_sec * 1000000 + stage2.tv_usec) - (stage1.tv_sec * 1000000 + stage1.tv_usec));
  params->stage2_time = (double) walltime / 1000000.0;
  walltime = ((stage3.tv_sec * 1000000 + stage3.tv_usec) - (stage2.tv_sec * 1000000 + stage2.tv_usec));
  params->stage3_time = (double) walltime / 1000000.0;
  walltime = ((stage4.tv_sec * 1000000 + stage4.tv_usec) - (stage3.tv_sec * 1000000 + stage3.tv_usec));
  params->stage4_time = (double) walltime / 1000000.0;
  walltime = ((stage4.tv_sec * 1000000 + stage4.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
  params->total_time = (double) walltime / 1000000.0;
  return(0);
}
