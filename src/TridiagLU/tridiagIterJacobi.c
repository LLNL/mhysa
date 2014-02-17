#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

int tridiagIterJacobi(double *a,double *b,double *c,double *x,
              int n,int ns,void *r,void *m)
{
  TridiagLU  *context = (TridiagLU*) r;
  int        iter,d,i,NT;
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

  /* check for zero along the diagonal */
  for (i=0; i<n; i++) {
    for (d=0; d<ns; d++) {
      if (b[i*ns+d]*b[i*ns+d] < context->atol*context->atol) {
        fprintf(stderr,"Error in tridiagIterJacobi(): Encountered zero on main diagonal!\n");
        return(1);
      }
    }
  }

  double *rhs = (double*) calloc (ns*n, sizeof(double));
  for (i=0; i<n; i++) {
    for (d=0; d<ns; d++) {
      rhs[i*ns+d] = x[i*ns+d]; /* save a copy of the rhs */
      x[i*ns+d]  /= b[i*ns+d]; /* initial guess          */
    }
  }

  double *recvbufL, *recvbufR, *sendbufL, *sendbufR;
  recvbufL = (double*) calloc (ns, sizeof(double));
  recvbufR = (double*) calloc (ns, sizeof(double));
  sendbufL = (double*) calloc (ns, sizeof(double));
  sendbufR = (double*) calloc (ns, sizeof(double));

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
    for (d=0; d<ns; d++)  recvbufL[d] = recvbufR[d] = 0;
#ifndef serial
    MPI_Request req[4] =  {MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if (rank)             MPI_Irecv(recvbufL,ns,MPI_DOUBLE,rank-1,2,*comm,&req[0]);
    if (rank != nproc-1)  MPI_Irecv(recvbufR,ns,MPI_DOUBLE,rank+1,3,*comm,&req[1]);
    for (d=0; d<ns; d++)  { sendbufL[d] = x[d]; sendbufR[d] = x[(n-1)*ns+d]; }
    if (rank)             MPI_Isend(sendbufL,ns,MPI_DOUBLE,rank-1,3,*comm,&req[2]);
    if (rank != nproc-1)  MPI_Isend(sendbufR,ns,MPI_DOUBLE,rank+1,2,*comm,&req[3]);
#endif

    /* calculate error norm - interior */
    if (context->evaluate_norm) {
      norm = 0;
      for (i=1; i<n-1; i++) {
        for (d=0; d<ns; d++) {
          norm  += ( (a[i*ns+d]*x[(i-1)*ns+d] + b[i*ns+d]*x[i*ns+d] + c[i*ns+d]*x[(i+1)*ns+d] - rhs[i*ns+d])
                   * (a[i*ns+d]*x[(i-1)*ns+d] + b[i*ns+d]*x[i*ns+d] + c[i*ns+d]*x[(i+1)*ns+d] - rhs[i*ns+d]) );
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
          norm  += ( (a[d]*recvbufL[d] + b[d]*x[d] + c[d]*x[d+ns*1]- rhs[d])
                   * (a[d]*recvbufL[d] + b[d]*x[d] + c[d]*x[d+ns*1]- rhs[d]) );
        }
        for (d=0; d<ns; d++) {
          norm  += ( (a[d+ns*(n-1)]*x[d+ns*(n-2)] + b[d+ns*(n-1)]*x[d+ns*(n-1)] + c[d+ns*(n-1)]*recvbufR[d] - rhs[d+ns*(n-1)])
                   * (a[d+ns*(n-1)]*x[d+ns*(n-2)] + b[d+ns*(n-1)]*x[d+ns*(n-1)] + c[d+ns*(n-1)]*recvbufR[d] - rhs[d+ns*(n-1)]) );
        }
      } else {
        for (d=0; d<ns; d++) {
          norm  += ( (a[d]*recvbufL[d] + b[d]*x[d] + c[d]*recvbufR[d] - rhs[d])
                   * (a[d]*recvbufL[d] + b[d]*x[d] + c[d]*recvbufR[d] - rhs[d]) );
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
        i = 0;    x[i*ns+d] = (rhs[i*ns+d] - a[i*ns+d]*recvbufL[d] - c[i*ns+d]*x[d+ns*(i+1)]  ) / b[i*ns+d];
        i = n-1;  x[i*ns+d] = (rhs[i*ns+d] - a[i*ns+d]*x[d+ns*(i-1)]   - c[i*ns+d]*recvbufR[d]) / b[i*ns+d];
      }
      for (i=1; i<n-1; i++) {
        for (d=0; d<ns; d++) {
          x[i*ns+d] = (rhs[i*ns+d] - a[i*ns+d]*x[d+ns*(i-1)] - c[i*ns+d]*x[d+ns*(i+1)]) / b[i*ns+d];
        }
      }
    } else {
      for (d=0; d<ns; d++) x[d] = (rhs[d] - a[d]*recvbufL[d] - c[d]*recvbufR[d]) / b[d];
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
