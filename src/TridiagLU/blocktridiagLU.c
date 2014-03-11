#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>
#include <matops.h>

int blocktridiagLU(double *a,double *b,double *c,double *x,
                   int n,int ns,int bs,void *r,void *m)
{
  TridiagLU       *params = (TridiagLU*) r;
  int             d,i,j,istart,iend,size;
  int             rank,nproc,bs2=bs*bs,nsbs=ns*bs;
  struct timeval  start,stage1,stage2,stage3,stage4;

#ifdef serial
  rank  = 0;
  nproc = 1;
#else
  MPI_Comm        *comm = (MPI_Comm*) m;
  const int       nvar = 4;
  int             ierr = 0;

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

  if ((ns == 0) || (n == 0) || (bs == 0)) return(0);
  double *xs1, *xp1;
  xs1 = (double*) calloc (nsbs, sizeof(double));
  xp1 = (double*) calloc (nsbs, sizeof(double));
  for (i=0; i<nsbs; i++) xs1[i] = xp1[i] = 0;

  /* Stage 1 - Parallel elimination of subdiagonal entries */
  istart  = (rank == 0 ? 1 : 2);
  iend    = n;
  for (i = istart; i < iend; i++) {
    double binv[bs2], factor[bs2];
    for (d = 0; d < ns; d++) {
      _MatrixInvert_           (b+((i-1)*ns+d)*bs2,binv,bs);
      _MatrixMultiply_         (a+(i*ns+d)*bs2,binv,factor,bs);
      _MatrixMultiplySubtract_ (b+(i*ns+d)*bs2,factor,c+((i-1)*ns+d)*bs2,bs);
      _MatrixZero_             (a+(i*ns+d)*bs2,bs);
      _MatrixMultiplySubtract_ (a+(i*ns+d)*bs2,factor,a+((i-1)*ns+d)*bs2,bs);
      _MatVecMultiplySubtract_ (x+(i*ns+d)*bs ,factor,x+((i-1)*ns+d)*bs ,bs);
      if (rank) {
        _MatrixMultiply_         (c+d*bs2,binv,factor,bs);
        _MatrixZero_             (c+d*bs2,bs);
        _MatrixMultiplySubtract_ (c+d*bs2,factor,c+((i-1)*ns+d)*bs2,bs);
        _MatrixMultiplySubtract_ (b+d*bs2,factor,a+((i-1)*ns+d)*bs2,bs);
        _MatVecMultiplySubtract_ (x+d*bs ,factor,x+((i-1)*ns+d)*bs ,bs);
      }
    }
  }

  /* end of stage 1 */
  gettimeofday(&stage1,NULL);

  /* Stage 2 - Eliminate the first sub- & super-diagonal entries */
  /* This needs the last (a,b,c,x) from the previous process     */
#ifndef serial
  double *sendbuf, *recvbuf;
  size = ns*bs2*(nvar-1)+nsbs;
  sendbuf = (double*) calloc (size, sizeof(double));
  recvbuf = (double*) calloc (size, sizeof(double));
  for (d=0; d<ns; d++) {
    for (i=0; i<bs2; i++) {
      sendbuf[(d*nvar+0)*bs2+i] = a[((n-1)*ns+d)*bs2+i]; 
      sendbuf[(d*nvar+1)*bs2+i] = b[((n-1)*ns+d)*bs2+i]; 
      sendbuf[(d*nvar+2)*bs2+i] = c[((n-1)*ns+d)*bs2+i];
    }
    for (i=0; i<bs; i++) sendbuf[(d*nvar+3)*bs2+i] = x[((n-1)*ns+d)*bs+i];
  }
  if (nproc > 1) {
    MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if (rank)             MPI_Irecv(recvbuf,size,MPI_DOUBLE,rank-1,1436,*comm,&req[0]);
    if (rank != nproc-1)  MPI_Isend(sendbuf,size,MPI_DOUBLE,rank+1,1436,*comm,&req[1]);
    MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  }
  /* The first process sits this one out */
  if (rank) {
    for (d = 0; d < ns; d++) {
      double am1[bs2], bm1[bs2], cm1[bs2], xm1[bs];
      for (i=0; i<bs2; i++) {
        am1[i] = recvbuf[(d*nvar+0)*bs2+i]; 
        bm1[i] = recvbuf[(d*nvar+1)*bs2+i]; 
        cm1[i] = recvbuf[(d*nvar+2)*bs2+i]; 
      }
      for (i=0; i<bs; i++) xm1[i] = recvbuf[(d*nvar+3)*bs2+i];
      double factor[bs2], binv[bs2];
      _MatrixInvert_           (bm1,binv,bs);
      _MatrixMultiply_         (a+d*bs2,binv,factor,bs);
      _MatrixMultiplySubtract_ (b+d*bs2,factor,cm1,bs);
      _MatrixZero_             (a+d*bs2,bs);
      _MatrixMultiplySubtract_ (a+d*bs2,factor,am1,bs);
      _MatVecMultiplySubtract_ (x+d*bs ,factor,xm1,bs);
      
      _MatrixInvert_           (b+((n-1)*ns+d)*bs2,binv,bs); if (ierr) return(ierr);
      _MatrixMultiply_         (c+d*bs2,binv,factor,bs);
      _MatrixMultiplySubtract_ (b+d*bs2,factor,a+((n-1)*ns+d)*bs2,bs);
      _MatrixZero_             (c+d*bs2,bs);
      _MatrixMultiplySubtract_ (c+d*bs2,factor,c+((n-1)*ns+d)*bs2,bs);
      _MatVecMultiplySubtract_ (x+d*bs ,factor,x+((n-1)*ns+d)*bs ,bs);
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
    double *zero, *eye;
    zero = (double*) calloc (ns*bs2, sizeof(double));
    eye  = (double*) calloc (ns*bs2, sizeof(double));
    for (d=0; d<ns*bs2; d++) zero[d] = eye[d] = 0.0;
    for (d=0; d<ns; d++) {
      for (i=0; i<bs; i++) eye[d*bs2+(i*bs+i)] = 1.0;
    }

    if (!strcmp(params->reducedsolvetype,_TRIDIAG_GS_)) {
      /* not supported */
      fprintf(stderr,"Error in blocktridiagLU(): Gather-and-solve for reduced system not available.\n");
      return(1);
    } else if (!strcmp(params->reducedsolvetype,_TRIDIAG_JACOBI_)) {
      /* Solving the reduced system iteratively with the Jacobi method */
      if (rank) ierr = blocktridiagIterJacobi(a,b,c,x,1,ns,bs,params,comm);
      else      ierr = blocktridiagIterJacobi(zero,eye,zero,zero,1,ns,bs,params,comm);
    }
    free(zero);
    free(eye);

    /* Each process, get the first x of the next process */
    MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    for (d=0; d<nsbs; d++)  xs1[d] = x[d];
    if (rank+1 < nproc) MPI_Irecv(xp1,nsbs,MPI_DOUBLE,rank+1,1323,*comm,&req[0]);
    if (rank)           MPI_Isend(xs1,nsbs,MPI_DOUBLE,rank-1,1323,*comm,&req[1]);
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
    double binv[bs2],xt[bs];
    _MatrixInvert_           (b+(istart*ns+d)*bs2,binv,bs);
    _MatVecMultiplySubtract_ (x+(istart*ns+d)*bs ,a+(istart*ns+d)*bs2,x  +d*bs,bs);
    _MatVecMultiplySubtract_ (x+(istart*ns+d)*bs ,c+(istart*ns+d)*bs2,xp1+d*bs,bs);
    _MatVecMultiply_         (binv,x+(istart*ns+d)*bs,xt,bs);
    for (j=0; j<bs; j++)   x[(istart*ns+d)*bs+j]=xt[j];
  }
  for (i = istart-1; i > iend-1; i--) {
    for (d = 0; d < ns; d++) {
      double binv[bs2],xt[bs];
      _MatrixInvert_           (b+(i*ns+d)*bs2,binv,bs);
      _MatVecMultiplySubtract_ (x+(i*ns+d)*bs ,c+(i*ns+d)*bs2,x+((i+1)*ns+d)*bs,bs);
      _MatVecMultiplySubtract_ (x+(i*ns+d)*bs ,a+(i*ns+d)*bs2,x+d*bs,bs);
      _MatVecMultiply_         (binv,x+(i*ns+d)*bs,xt,bs);
      for (j=0; j<bs; j++)   x[(i*ns+d)*bs+j] = xt[j];
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
