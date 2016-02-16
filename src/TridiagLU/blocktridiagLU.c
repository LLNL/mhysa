/*! @file blocktridiagLU.c
    @brief Solve block tridiagonal systems using the LU decomposition method
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>
#include <matops.h>

/*!
  Solve block tridiagonal (non-periodic) systems of equations using parallel LU decomposition: 
  This function can solve multiple independent systems with one call. The systems need not share 
  the same left- or right-hand-sides. The iterative substructuring method is used in this
  function that can be briefly described through the following 4 stages:
  + Stage 1: Parallel elimination of the tridiagonal blocks on each processor comprising all 
    points of the subdomain except the 1st point (unless its the 1st global point, i.e., a 
    physical boundary)
  + Stage 2: Elimination of the 1st row on each processor (except the 1st processor) using the 
    last row of the previous processor.
  + Stage 3: Solution of the reduced tridiagonal system that represents the coupling of the
    system across the processors, using blocktridiagIterJacobi() in this implementation.
  + Stage 4: Backward-solve to obtain the final solution

  Specific details of the method implemented here are available in:
  + Ghosh, D., Constantinescu, E. M., Brown, J., "Scalable Nonlinear Compact Schemes", 
    Technical Memorandum, ANL/MCS-TM-340, Argonne National Laboratory, April 2014,
    (http://www.mcs.anl.gov/publication/scalable-nonlinear-compact-schemes)
    (also available at http://debog.github.io/Files/2014_Ghosh_Consta_Brown_MCSTR340.pdf).
  + Ghosh, D., Constantinescu, E. M., Brown, J., Efficient Implementation of Nonlinear 
    Compact Schemes on Massively Parallel Platforms, SIAM Journal on Scientific Computing, 
    37 (3), 2015, C354–C383 (http://dx.doi.org/10.1137/140989261).

  More references on this class of parallel tridiagonal solvers:
  + E. Polizzi and A. H. Sameh, "A parallel hybrid banded system solver: The SPIKE algorithm",
    Parallel Comput., 32 (2006), pp. 177–194.
  + E. Polizzi and A. H. Sameh, "SPIKE: A parallel environment for solving banded linear systems",
    Comput. & Fluids, 36 (2007), pp. 113–120.

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
int blocktridiagLU(
                    double  *a, /*!< Array containing the sub-diagonal elements */
                    double  *b, /*!< Array containing the diagonal elements */
                    double  *c, /*!< Array containing the super-diagonal elements */
                    double  *x, /*!< Right-hand side; will contain the solution on exit */
                    int     n,  /*!< Local size of the system on this processor (*not*
                                     multiplied by the block size) */
                    int     ns, /*!< Number of systems to solve */
                    int     bs, /*!< Block size */
                    void    *r, /*!< Object of type #TridiagLU (contains wall times at exit) */
                    void    *m  /*!< MPI communicator */
                  )
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
      sendbuf[(0*ns+d)*bs2+i] = a[((n-1)*ns+d)*bs2+i]; 
      sendbuf[(1*ns+d)*bs2+i] = b[((n-1)*ns+d)*bs2+i]; 
      sendbuf[(2*ns+d)*bs2+i] = c[((n-1)*ns+d)*bs2+i];
    }
    for (i=0; i<bs; i++) sendbuf[3*ns*bs2+d*bs+i] = x[((n-1)*ns+d)*bs+i];
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
        am1[i] = recvbuf[(0*ns+d)*bs2+i]; 
        bm1[i] = recvbuf[(1*ns+d)*bs2+i]; 
        cm1[i] = recvbuf[(2*ns+d)*bs2+i]; 
      }
      for (i=0; i<bs; i++) xm1[i] = recvbuf[3*ns*bs2+d*bs+i];
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
