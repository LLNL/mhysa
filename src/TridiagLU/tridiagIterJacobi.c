/*! @file tridiagIterJacobi.c
    @brief Solve tridiagonal systems of equations with the point Jacobi method
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

/*!
  Solve tridiagonal (non-periodic) systems of equations using point Jacobi iterations: 
  This function can solve multiple independent systems with one call. The systems need not share 
  the same left- or right-hand-sides. The initial guess is taken as the solution of 
  \f{equation}{
    {\rm diag}\left[{\bf b}\right]{\bf x} = {\bf r}
  \f}
  where \f${\bf b}\f$ represents the diagonal elements of the tridiagonal system, and 
  \f${\bf r}\f$ is the right-hand-side, stored in \f${\bf x}\f$ at the start of this
  function.

  Array layout: The arguments \a a, \a b, \a c, and \a x are local 1D arrays (containing
  this processor's part of the subdiagonal, diagonal, superdiagonal, and right-hand-side)
  of size (\a n X \a ns), where \a n is the local size of the system, and \a ns is
  the number of independent systems to solve. The ordering of the elements in these arrays 
  is as follows:
  + Elements of the same row for each of the independent systems are stored adjacent to each 
    other.

  For example, consider the following systems:
  \f{equation}{
    \left[\begin{array}{ccccc}
      b_0^k & c_0^k &               &       &       \\
      a_1^k & b_1^k & c_1^k         &       &       \\
            & a_2^k & b_2^k & c_2^k &       &       \\
            &       & a_3^k & b_3^k & c_3^k &       \\
            &       &       & a_4^k & b_4^k & c_4^k \\
    \end{array}\right]
    \left[\begin{array}{c} x_0^k \\ x_1^k \\ x_2^k \\ x_3^k \\ x_4^k \end{array}\right]
    =
    \left[\begin{array}{c} r_0^k \\ r_1^k \\ r_2^k \\ r_3^k \\ r_4^k \end{array}\right];
    \ \ k= 1,\cdots,ns
  \f}
  and let \f$ ns = 3\f$. Note that in the code, \f$x\f$ and \f$r\f$ are the same array \a x.
  
  Then, the array \a b must be a 1D array with the following layout of elements:\n
  [\n
  b_0^0, b_0^1, b_0^2, (diagonal element of the first row in each system) \n
  b_1^0, b_1^1, b_1^2, (diagonal element of the second row in each system) \n
  ..., \n
  b_{n-1}^0, b_{n-1}^1, b_{n-1}^2 (diagonal element of the last row in each system) \n
  ]\n
  The arrays \a a, \a c, and \a x are stored similarly. 
  
  Notes:
  + This function does *not* preserve the sub-diagonal, diagonal, super-diagonal elements
    and the right-hand-sides. 
  + The input array \a x contains the right-hand-side on entering the function, and the 
    solution on exiting it.
*/
int tridiagIterJacobi(
                        double  *a, /*!< Array containing the sub-diagonal elements */
                        double  *b, /*!< Array containing the diagonal elements */
                        double  *c, /*!< Array containing the super-diagonal elements */
                        double  *x, /*!< Right-hand side; will contain the solution on exit */
                        int     n,  /*!< Local size of the system on this processor */
                        int     ns, /*!< Number of systems to solve */
                        void    *r, /*!< Object of type #TridiagLU */
                        void    *m  /*!< MPI communicator */
                     )
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
