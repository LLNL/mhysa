#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>
#include <matops.h>

int blocktridiagIterJacobi(double *a,double *b,double *c,double *x,
                           int n,int ns,int bs,void *r,void *m)
{
  TridiagLU  *context = (TridiagLU*) r;
  int        iter,d,i,j,NT,bs2=bs*bs,nsbs=ns*bs;
  double     norm=0,norm0=0,global_norm=0;

#ifndef serial
  MPI_Comm  *comm = (MPI_Comm*) m;
  int       rank,nproc;

  if (comm) {
    MPI_Comm_size(*comm,&nproc);
    MPI_Comm_rank(*comm,&rank);
  } else {
    rank  = 0;
    nproc = 1;
  }
#endif

  if (!context) {
    fprintf(stderr,"Error in tridiagIterJacobi(): NULL pointer passed for parameters!\n");
    return(-1);
  }

  double *rhs = (double*) calloc (ns*n*bs, sizeof(double));
  for (i=0; i<ns*n*bs; i++) rhs[i] = x[i]; /* save a copy of the rhs */
  /* initial guess */
  for (i=0; i<n; i++) {
    for (d=0; d<ns; d++) {
      double binv[bs2];
      _MatrixInvert_    (b+(i*ns+d)*bs2,binv,bs);
      _MatVecMultiply_  (binv,rhs+(i*ns+d)*bs,x+(i*ns+d)*bs,bs);
    }
  }

  double *recvbufL, *recvbufR, *sendbufL, *sendbufR;
  recvbufL = (double*) calloc (nsbs, sizeof(double));
  recvbufR = (double*) calloc (nsbs, sizeof(double));
  sendbufL = (double*) calloc (nsbs, sizeof(double));
  sendbufR = (double*) calloc (nsbs, sizeof(double));

  /* total number of points */
#ifdef serial
  if (context->evaluate_norm)    NT = n;
  else                           NT = 0;
#else
  if (context->evaluate_norm) MPI_Allreduce(&n,&NT,1,MPI_INT,MPI_SUM,*comm);
  else NT = 0;
#endif

#ifdef serial
    if (context->verbose) printf("\n");
#else
    if (context->verbose && (!rank)) printf("\n");
#endif

  iter = 0;
  while(1) {

    /* evaluate break conditions */
    if (    (iter >= context->maxiter) 
        ||  (iter && context->evaluate_norm && (global_norm < context->atol)) 
        ||  (iter && context->evaluate_norm && (global_norm/norm0 < context->rtol))  ) {
      break;
    }

    /* Communicate the boundary x values between processors */
    for (d=0; d<nsbs; d++)  recvbufL[d] = recvbufR[d] = 0;
#ifndef serial
    MPI_Request req[4] =  {MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if (rank)             MPI_Irecv(recvbufL,nsbs,MPI_DOUBLE,rank-1,2,*comm,&req[0]);
    if (rank != nproc-1)  MPI_Irecv(recvbufR,nsbs,MPI_DOUBLE,rank+1,3,*comm,&req[1]);
    for (d=0; d<nsbs; d++) {
      sendbufL[d] = x[d];
      sendbufR[d] = x[(n-1)*nsbs+d];
    }
    if (rank)             MPI_Isend(sendbufL,nsbs,MPI_DOUBLE,rank-1,3,*comm,&req[2]);
    if (rank != nproc-1)  MPI_Isend(sendbufR,nsbs,MPI_DOUBLE,rank+1,2,*comm,&req[3]);
#endif

    /* calculate error norm - interior */
    if (context->evaluate_norm) {
      norm = 0;
      for (i=1; i<n-1; i++) {
        for (d=0; d<ns; d++) {
          double err[bs]; for (j=0; j<bs; j++) err[j] = rhs[(i*ns+d)*bs+j];
          _MatVecMultiplySubtract_(err,a+(i*ns+d)*bs2,x+((i-1)*ns+d)*bs,bs);
          _MatVecMultiplySubtract_(err,b+(i*ns+d)*bs2,x+((i  )*ns+d)*bs,bs);
          _MatVecMultiplySubtract_(err,c+(i*ns+d)*bs2,x+((i+1)*ns+d)*bs,bs);
          for (j=0; j<bs; j++) norm += (err[j]*err[j]);
        }
      }
    }
    /* calculate error norm - boundary */
#ifndef serial
    MPI_Waitall(4,req,MPI_STATUS_IGNORE);
#endif
    if (context->evaluate_norm) {
      if (n > 1) {
        for (d=0; d<ns; d++) {
          double err[bs]; for (j=0; j<bs; j++) err[j] = rhs[d*bs+j];
          _MatVecMultiplySubtract_(err,a+d*bs2,recvbufL+d*bs,bs);
          _MatVecMultiplySubtract_(err,b+d*bs2,x+d*bs,bs);
          _MatVecMultiplySubtract_(err,c+d*bs2,x+(ns+d)*bs,bs);
          for (j=0; j<bs; j++) norm += (err[j]*err[j]);
        }
        for (d=0; d<ns; d++) {
          double err[bs]; for (j=0; j<bs; j++) err[j] = rhs[(d+ns*(n-1))*bs+j];
          _MatVecMultiplySubtract_(err,a+(d+ns*(n-1))*bs2,x+(d+ns*(n-2))*bs,bs);
          _MatVecMultiplySubtract_(err,b+(d+ns*(n-1))*bs2,x+(d+ns*(n-1))*bs,bs);
          _MatVecMultiplySubtract_(err,c+(d+ns*(n-1))*bs2,recvbufR+d*bs,bs);
          for (j=0; j<bs; j++) norm += (err[j]*err[j]);
        }
      } else {
        for (d=0; d<ns; d++) {
          double err[bs]; for (j=0; j<bs; j++) err[j] = rhs[d*bs+j];
          _MatVecMultiplySubtract_(err,a+d*bs2,recvbufL+d*bs,bs);
          _MatVecMultiplySubtract_(err,b+d*bs2,x+d*bs,bs);
          _MatVecMultiplySubtract_(err,c+d*bs2,recvbufR+d*bs,bs);
          for (j=0; j<bs; j++) norm += (err[j]*err[j]);
        }
      }
      /* sum over all processes */
#ifdef serial
      global_norm = norm;
#else
      MPI_Allreduce(&norm,&global_norm,1,MPI_DOUBLE,MPI_SUM,*comm);
#endif
      global_norm = sqrt(global_norm/NT);
      if (!iter) norm0 = global_norm;
    } else {
      norm = -1.0;
      global_norm = -1.0;
    }

#ifdef serial
    if (context->verbose)
#else
    if (context->verbose && (!rank))
#endif
      printf("\t\titer: %d, norm: %1.16E\n",iter,global_norm);

    /* correct the solution for this iteration */
    if (n > 1) {
      for (d=0; d<ns; d++) {
        double xt[bs],binv[bs2];
        
        i = 0;    
        for (j=0; j<bs; j++) xt[j] = rhs[(i*ns+d)*bs+j];
        _MatVecMultiplySubtract_(xt,a+(i*ns+d)*bs2,recvbufL+d*bs,bs);
        _MatVecMultiplySubtract_(xt,c+(i*ns+d)*bs2,x+(d+ns*(i+1))*bs,bs);
        _MatrixInvert_(b+(i*ns+d)*bs2,binv,bs);
        _MatVecMultiply_(binv,xt,x+(i*ns+d)*bs,bs);

        i = n-1;  
        for (j=0; j<bs; j++) xt[j] = rhs[(i*ns+d)*bs+j];
        _MatVecMultiplySubtract_(xt,a+(i*ns+d)*bs2,x+(d+ns*(i-1))*bs,bs);
        _MatVecMultiplySubtract_(xt,c+(i*ns+d)*bs2,recvbufR+d*bs,bs);
        _MatrixInvert_(b+(i*ns+d)*bs2,binv,bs);
        _MatVecMultiply_(binv,xt,x+(i*ns+d)*bs,bs);
      }
      for (i=1; i<n-1; i++) {
        for (d=0; d<ns; d++) {
          double xt[bs],binv[bs2];
          for (j=0; j<bs; j++) xt[j] = rhs[(i*ns+d)*bs+j];
          _MatVecMultiplySubtract_(xt,a+(i*ns+d)*bs2,x+(d+ns*(i-1))*bs,bs);
          _MatVecMultiplySubtract_(xt,c+(i*ns+d)*bs2,x+(d+ns*(i+1))*bs,bs);
          _MatrixInvert_(b+(i*ns+d)*bs2,binv,bs);
          _MatVecMultiply_(binv,xt,x+(i*ns+d)*bs,bs);
        }
      }
    } else {
      for (d=0; d<ns; d++) {
        double xt[bs],binv[bs2];
        for (j=0; j<bs; j++) xt[j] = rhs[d*bs+j];
        _MatVecMultiplySubtract_(xt,a+d*bs2,recvbufL+d*bs,bs);
        _MatVecMultiplySubtract_(xt,c+d*bs2,recvbufR+d*bs,bs);
        _MatrixInvert_(b+d*bs2,binv,bs);
        _MatVecMultiply_(binv,xt,x+d*bs,bs);
      }
    }

    /* finished with this iteration */
    iter++;
  }

  /* save convergence information */
  context->exitnorm = (context->evaluate_norm ? global_norm : -1.0);
  context->exititer = iter;

  free(rhs);
  free(sendbufL);
  free(sendbufR);
  free(recvbufL);
  free(recvbufR);

  return(0);
}
