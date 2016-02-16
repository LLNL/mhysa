/*! @file blocktridiagIterJacobi.c
    @brief Solve a block tridiagonal system with the Jacobi method
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>
#include <matops.h>

/*!
  Solve block tridiagonal (non-periodic) systems of equations using point Jacobi iterations: 
  This function can solve multiple independent systems with one call. The systems need not share 
  the same left- or right-hand-sides. The initial guess is taken as the solution of 
  \f{equation}{
    {\rm diag}\left[{\bf b}\right]{\bf x} = {\bf r}
  \f}
  where \f${\bf b}\f$ represents the diagonal elements of the tridiagonal system, and 
  \f${\bf r}\f$ is the right-hand-side, stored in \f${\bf x}\f$ at the start of this
  function.

  Array layout: The arguments \a a, \a b, and \a c are local 1D arrays (containing
  this processor's part of the subdiagonal, diagonal, and superdiagonal)
  of size (\a n X \a ns X \a bs^2), and \a x is a local 1D array (containing this
  processor's part of the right-hand-side, and will contain the solution on exit) 
  of size (\a n X \a ns X \a bs), where \a n is the local size of the system, \a ns is
  the number of independent systems to solve, and \a bs is the block size. The ordering
  of the elements in these arrays is as follows:
  + Each block is stored in the row-major format.
  + Blocks of the same row for each of the independent systems are stored adjacent to each 
    other.

  For example, consider the following systems:
  \f{equation}{
    \left[\begin{array}{ccccc}
      B_0^k & C_0^k &               &       &       \\
      A_1^k & B_1^k & C_1^k         &       &       \\
            & A_2^k & B_2^k & C_2^k &       &       \\
            &       & A_3^k & B_3^k & C_3^k &       \\
            &       &       & A_4^k & B_4^k & C_4^k \\
    \end{array}\right]
    \left[\begin{array}{c} X_0^k \\ X_1^k \\ X_2^k \\ X_3^k \\ X_4^k \end{array}\right]
    =
    \left[\begin{array}{c} R_0^k \\ R_1^k \\ R_2^k \\ R_3^k \\ R_4^k \end{array}\right];
    \ \ k= 1,\cdots,ns
  \f}
  where \f$A\f$, \f$B\f$, and \f$C\f$ are matrices of size \a bs = 2 (say), and let
  \f$ ns = 3\f$. In the equation above, we have
  \f{equation}{
    B_i^k = \left[\begin{array}{cc} b_{00,i}^k & b_{01,i}^k \\ b_{10,i}^k & b_{11,i}^k \end{array}\right],
    X_i^k = \left[\begin{array}{c} x_{0,i}^k \\ x_{1,i}^k \end{array} \right],
    R_i^k = \left[\begin{array}{c} r_{0,i}^k \\ r_{1,i}^k \end{array} \right]
  \f}
  Note that in the code, \f$X\f$ and \f$R\f$ are the same array \a x.
  
  Then, the array \a b must be a 1D array with the following layout of elements:\n
  [\n
  b_{00,0}^0, b_{01,0}^0, b_{10,0}^0, b_{11,0}^0, b_{00,0}^1, b_{01,0}^1, b_{10,0}^1, b_{11,0}^1,
  b_{00,0}^2, b_{01,0}^2, b_{10,0}^2, b_{11,0}^2, \n
  b_{00,1}^0, b_{01,1}^0, b_{10,1}^0, b_{11,1}^0, b_{00,1}^1, b_{01,1}^1, b_{10,1}^1, b_{11,1}^1, 
  b_{00,1}^2, b_{01,1}^2, b_{10,1}^2, b_{11,1}^2, \n
  ..., \n
  b_{00,n-1}^0, b_{01,n-1}^0, b_{10,n-1}^0, b_{11,n-1}^0, b_{00,n-1}^1, b_{01,n-1}^1, b_{10,n-1}^1, b_{11,n-1}^1, 
  b_{00,n-1}^2, b_{01,n-1}^2, b_{10,n-1}^2, b_{11,n-1}^2\n
  ]\n
  The arrays \a a and \a c are stored similarly. 
  
  The array corresponding to a vector (the solution and the right-hand-side \a x) must be a 1D array with the following 
  layout of elements:\n
  [\n
  x_{0,0}^0, x_{1,0}^0, x_{0,0}^1, x_{1,0}^1,x_{0,0}^2, x_{1,0}^2,\n
  x_{0,1}^0, x_{1,1}^0, x_{0,1}^1, x_{1,1}^1,x_{0,1}^2, x_{1,1}^2,\n
  ..., \n
  x_{0,n-1}^0, x_{1,n-1}^0, x_{0,n-1}^1, x_{1,n-1}^1,x_{0,n-1}^2, x_{1,n-1}^2\n
  ]\n


  Notes:
  + This function does *not* preserve the sub-diagonal, diagonal, super-diagonal elements
    and the right-hand-sides. 
  + The input array \a x contains the right-hand-side on entering the function, and the 
    solution on exiting it.
*/
int blocktridiagIterJacobi(
                            double  *a, /*!< Array containing the sub-diagonal elements */
                            double  *b, /*!< Array containing the diagonal elements */
                            double  *c, /*!< Array containing the super-diagonal elements */
                            double  *x, /*!< Right-hand side; will contain the solution on exit */
                            int     n,  /*!< Local size of the system on this processor (*not*
                                             multiplied by the block size) */
                            int     ns, /*!< Number of systems to solve */
                            int     bs, /*!< Block size */
                            void    *r, /*!< Object of type #TridiagLU */
                            void    *m  /*!< MPI communicator */
                          )
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
