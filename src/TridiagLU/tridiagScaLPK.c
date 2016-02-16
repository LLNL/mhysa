/*! @file tridiagScaLPK.c
    @brief Wrapper function to solve tridiagonal systems using ScaLAPACK's pddtsv
    @author Debojyoti Ghosh
*/

#ifdef with_scalapack

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

/*! 
  ScaLAPACK's parallel tridiagonal solver: See 
  http://www.netlib.org/scalapack/explore-html/d7/d86/pddtsv_8f_source.html
  for its documentation.
*/
extern void pddtsv_();

/*!
  Solve tridiagonal (non-periodic) systems of equations using ScaLAPACK's pddtsv:
  This function can solve multiple independent systems with one call. The systems need not share 
  the same left- or right-hand-sides. 
  + This function is compiled only with the compilation flag "-Dwith_scalapack" is specified.
  + This function calls the ScaLAPACK function for solving tridiagonal systems individually
    for each system, and thus may not be efficient.

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
int tridiagScaLPK(
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
  TridiagLU       *params = (TridiagLU*) r;
  int             rank,nproc,nglobal,nrhs,i,s,ia,ib,desca[9],descb[9],ierr,
                  lwork;
  double          *dl,*d,*du,*rhs,*work;
  struct timeval  start,end;

#ifdef serial
  rank  = 0;
  nproc = 1;
  nglobal=n;
#else
  MPI_Comm        *comm = (MPI_Comm*) m;

  if (comm) {
    MPI_Comm_size(*comm,&nproc);
    MPI_Comm_rank(*comm,&rank);
  } else {
    rank  = 0;
    nproc = 1;
  }
  MPI_Allreduce(&n,&nglobal,1,MPI_INT,MPI_SUM,*comm);
#endif

  /* check */
  if (nglobal%n != 0) {
    if (!rank) {
      fprintf(stderr,"Error: The ScaLAPACK wrapper can only handle cases where the global ");
      fprintf(stderr,"size of system is an integer multiple of no. of processes.\n");
    }
    return(1);
  }

  if (!params) {
    fprintf(stderr,"Error in tridiagLU(): NULL pointer passed for parameters.\n");
    return(1);
  }


  nrhs = 1;
  ia = 1;
  ib = 1;

  lwork = (12*nproc+3*n) + ( (8*nproc) > (10*nproc+4*nrhs) ? (8*nproc) : (10*nproc+4*nrhs) );

  desca[0] = 501;
  desca[1] = params->blacs_ctxt;
  desca[2] = nglobal;
  desca[3] = n;
  desca[4] = 0;
  desca[5] = 0;
  desca[6] = 0;
  desca[7] = 0;
  desca[8] = 0;

  descb[0] = 502;
  descb[1] = params->blacs_ctxt;
  descb[2] = nglobal;
  descb[3] = n;
  descb[4] = 0;
  descb[5] = n;
  descb[6] = 0;
  descb[7] = 0;
  descb[8] = 0;

  dl    = (double*) calloc (n,sizeof(double));
  d     = (double*) calloc (n,sizeof(double));
  du    = (double*) calloc (n,sizeof(double));
  rhs   = (double*) calloc (n,sizeof(double));
  work  = (double*) calloc(lwork,sizeof(double));

  params->total_time = 0.0;
  params->stage1_time = 0.0;
  params->stage2_time = 0.0;
  params->stage3_time = 0.0;
  params->stage4_time = 0.0;
  for (s=0; s<ns; s++) {

    for (i=0; i<n; i++) {
      dl[i] = a[i*ns+s];
      d [i] = b[i*ns+s];
      du[i] = c[i*ns+s];
      rhs[i]= x[i*ns+s];
    }

    /* call the ScaLAPACK function */
    gettimeofday(&start,NULL);
    pddtsv_(&nglobal,&nrhs,dl,d,du,&ia,desca,rhs,&ib,descb,work,&lwork,&ierr);
    gettimeofday(&end,NULL);
    if (ierr) return(ierr);
    
    for (i=0; i<n; i++) x[i*ns+s] = rhs[i];

    long long walltime;
    walltime = ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
    params->total_time += (double) walltime / 1000000.0;
  }


  free(dl);
  free(d);
  free(du);
  free(rhs);
  free(work);

  return(0);
}

#endif
