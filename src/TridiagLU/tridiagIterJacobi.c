#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

int tridiagIterJacobi(double **a,double **b,double **c,double **x,
              int n,int ns,void *r,void *m)
{
  TridiagLU  *context = (TridiagLU*) r;
  int        iter,d,i,NT;
  double     norm=0,norm0=0,global_norm=0,**rhs;
  double     *sendbufL,*recvbufL;
  double     *sendbufR,*recvbufR;

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
    fprintf(stderr,"Error in tridiagIterJacobi(): NULL pointer passed for paramters!\n");
    return(-1);
  }

  /* check for zero along the diagonal */
  for (d=0; d<ns; d++) {
    for (i=0; i<n; i++) {
      if (b[d][i]*b[d][i] < context->atol*context->atol) {
        fprintf(stderr,"Error in tridiagIterJacobi(): Encountered zero on main diagonal!\n");
        return(1);
      }
    }
  }

  rhs = (double**) calloc (ns,sizeof(double*));
  for (d=0; d<ns; d++) {
    rhs[d] = (double*) calloc (n,sizeof(double));
    for (i=0; i<n; i++) {
      rhs[d][i] = x[d][i]; /* save a copy of the rhs */
      x[d][i]  /= b[d][i]; /* initial guess          */
    }
  }

  recvbufL = (double*) calloc (ns,sizeof(double));
  recvbufR = (double*) calloc (ns,sizeof(double));
  sendbufL = (double*) calloc (ns,sizeof(double));
  sendbufR = (double*) calloc (ns,sizeof(double));

  /* total number of points */
  if (context->evaluate_norm)
#ifdef serial
    NT = n;
#else
    MPI_Allreduce(&n,&NT,1,MPI_INT,MPI_SUM,*comm);
#endif

#ifdef serial
    if (context->verbose)
#else
    if (context->verbose && (!rank))
#endif
      printf("\n");

  iter = 0;
  while(1) {

    /* Communicate the boundary x values between processors */
    for (d=0; d<ns; d++)  recvbufL[d] = recvbufR[d] = 0;
#ifndef serial
    MPI_Request req[4] =  {MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    if (rank)             MPI_Irecv(recvbufL,ns,MPI_DOUBLE,rank-1,2,*comm,&req[0]);
    if (rank != nproc-1)  MPI_Irecv(recvbufR,ns,MPI_DOUBLE,rank+1,3,*comm,&req[1]);
    for (d=0; d<ns; d++)  { sendbufL[d] = x[d][0]; sendbufR[d] = x[d][n-1]; }
    if (rank)             MPI_Isend(sendbufL,ns,MPI_DOUBLE,rank-1,3,*comm,&req[2]);
    if (rank != nproc-1)  MPI_Isend(sendbufR,ns,MPI_DOUBLE,rank+1,2,*comm,&req[3]);
#endif

    /* calculate error norm - interior */
    if (context->evaluate_norm) {
      norm = 0;
      for (d=0; d<ns; d++) {
        for (i=1; i<n-1; i++) {
          norm  += ( (a[d][i]*x[d][i-1] + b[d][i]*x[d][i] + c[d][i]*x[d][i+1] - rhs[d][i])
                   * (a[d][i]*x[d][i-1] + b[d][i]*x[d][i] + c[d][i]*x[d][i+1] - rhs[d][i]) );
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
          norm  += ( (a[d][0]*recvbufL[d] + b[d][0]*x[d][0] + c[d][0]*x[d][1]- rhs[d][0])
                   * (a[d][0]*recvbufL[d] + b[d][0]*x[d][0] + c[d][0]*x[d][1]- rhs[d][0]) );
        }
        for (d=0; d<ns; d++) {
          norm  += ( (a[d][n-1]*x[d][n-2] + b[d][n-1]*x[d][n-1] + c[d][n-1]*recvbufR[d] - rhs[d][n-1])
                   * (a[d][n-1]*x[d][n-2] + b[d][n-1]*x[d][n-1] + c[d][n-1]*recvbufR[d] - rhs[d][n-1]) );
        }
      } else {
        for (d=0; d<ns; d++) {
          norm  += ( (a[d][0]*recvbufL[d] + b[d][0]*x[d][0] + c[d][0]*recvbufR[d] - rhs[d][0])
                   * (a[d][0]*recvbufL[d] + b[d][0]*x[d][0] + c[d][0]*recvbufR[d] - rhs[d][0]) );
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

    /* evaluate break conditions */
    if (    (iter > context->maxiter) 
        ||  (context->evaluate_norm && (global_norm < context->atol)) 
        ||  (context->evaluate_norm && (global_norm/norm0 < context->rtol))  ) {
      break;
    }

    /* correct the solution for this iteration */
    for (d=0; d<ns; d++) {
      if (n > 1) {
        i = 0;    x[d][i] = (rhs[d][i] - a[d][i]*recvbufL[d] - c[d][i]*x[d][i+1]  ) / b[d][i];
        i = n-1;  x[d][i] = (rhs[d][i] - a[d][i]*x[d][i-1]   - c[d][i]*recvbufR[d]) / b[d][i];
        for (i=1; i<n-1; i++) x[d][i] = (rhs[d][i] - a[d][i]*x[d][i-1] - c[d][i]*x[d][i+1]) / b[d][i];
      } else x[d][0] = (rhs[d][0] - a[d][0]*recvbufL[d] - c[d][0]*recvbufR[d]) / b[d][0];
    }

    /* finished with this iteration */
    iter++;
  }

  /* save convergence information */
  context->exitnorm = (context->evaluate_norm ? global_norm : -1.0);
  context->exititer = iter;

  free(recvbufL);
  free(recvbufR);
  free(sendbufL);
  free(sendbufR);
  for (d=0; d<ns; d++) free(rhs[d]); free(rhs);

  return(0);
}
