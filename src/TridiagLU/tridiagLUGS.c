/*! @file tridiagLUGS.c
    @brief Solve tridiagonal systems of equations in parallel using "gather-and-solve"
    @author Debojyoti Ghosh
*/

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

/*!
  Solve tridiagonal (non-periodic) systems of equations using the gather-and-solve approach: 
  This function can solve multiple independent systems with one call. The systems need not share 
  the same left- or right-hand-sides. The "gather-and-solve" approach gathers a tridiagonal
  system on one processor and solves it using tridiagLU() (sending NULL as the argument for 
  MPI communicator to indicate that a serial solution is desired). Given multiple tridiagonal
  systems (\a ns > 1), this function will gather the systems on different processors in an
  optimal way, and thus each processor will solve some of the systems. After the system(s) is
  (are) solved, the solution(s) is (are) scattered back to the original processors.

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
int tridiagLUGS(
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
  TridiagLU *context = (TridiagLU*) r;
  if (!context) {
    fprintf(stderr,"Error in tridiagLUGS(): NULL pointer passed for parameters.\n");
    return(-1);
  }
#ifdef serial

  /* Serial compilation */
  return(tridiagLU(a,b,c,x,n,ns,context,m));

#else

  int         d,i,ierr = 0,dstart,istart,p,q;
  const int   nvar = 4;
  double      *sendbuf,*recvbuf;
  int         rank,nproc;

  /* Parallel compilation */
  MPI_Comm  *comm = (MPI_Comm*) m;
  if (!comm) return(tridiagLU(a,b,c,x,n,ns,context,NULL));
  MPI_Comm_size(*comm,&nproc);
  MPI_Comm_rank(*comm,&rank);

  if ((ns == 0) || (n == 0)) return(0);

  /* 
    each process needs to know the local sizes of every other process 
    and total size of the system
  */
  int *N = (int*) calloc (nproc, sizeof(int));
  MPI_Allgather(&n,1,MPI_INT,N,1,MPI_INT,*comm);
  int NT = 0; for (i=0; i<nproc; i++) NT += N[i];

  /* counts and displacements for gather and scatter operations */
  int *counts = (int*) calloc (nproc, sizeof(int));
  int *displ  = (int*) calloc (nproc, sizeof(int));

  /* on all processes, calculate the number of systems each     */
  /* process has to solve                                       */
  int *ns_local = (int*) calloc (nproc,sizeof(int));
  for (p=0; p<nproc; p++)    ns_local[p] = ns / nproc; 
  for (p=0; p<ns%nproc; p++) ns_local[p]++;

  /* allocate the arrays for the gathered tridiagonal systems */
  double *ra=NULL,*rb=NULL,*rc=NULL,*rx=NULL; 
  if (ns_local[rank] > 0) {
    ra = (double*) calloc (ns_local[rank]*NT,sizeof(double));
    rb = (double*) calloc (ns_local[rank]*NT,sizeof(double));
    rc = (double*) calloc (ns_local[rank]*NT,sizeof(double));
    rx = (double*) calloc (ns_local[rank]*NT,sizeof(double));
  }

  /* Gather the systems on each process */
  /* allocate receive buffer */
  if (ns_local[rank] > 0) 
    recvbuf = (double*) calloc (ns_local[rank]*nvar*NT,sizeof(double));
  else recvbuf = NULL;
  dstart = 0;
  for (p = 0; p < nproc; p++) {
    if (ns_local[p] > 0) {
      /* allocate send buffer and form the send packet of data */
      sendbuf = (double*) calloc (nvar*n*ns_local[p],sizeof(double));
      for (d = 0; d < ns_local[p]; d++) {
        for (i = 0; i < n; i++) {
          sendbuf[n*nvar*d+n*0+i] = a[d+dstart+ns*i];
          sendbuf[n*nvar*d+n*1+i] = b[d+dstart+ns*i];
          sendbuf[n*nvar*d+n*2+i] = c[d+dstart+ns*i];
          sendbuf[n*nvar*d+n*3+i] = x[d+dstart+ns*i];
        }
      }
      dstart += ns_local[p];

      /* gather these reduced systems on process with rank = p */
      if (rank == p) {
        for (q = 0; q < nproc; q++) {
          counts[q] = nvar*N[q]*ns_local[p];
          displ [q] = (q == 0 ? 0 : displ[q-1]+counts[q-1]);
        }
      }
      MPI_Gatherv(sendbuf,nvar*n*ns_local[p],MPI_DOUBLE,
                  recvbuf,counts,displ,MPI_DOUBLE,p,*comm);

      /* deallocate send buffer */
      free(sendbuf);
    }
  }
  /* extract the data from the recvbuf and solve */
  istart = 0;
  for (q = 0; q < nproc; q++) {
    for (d = 0; d < ns_local[rank]; d++) {
      for (i = 0; i < N[q]; i++) {
        ra[d+ns_local[rank]*(istart+i)] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+0*N[q]+i];
        rb[d+ns_local[rank]*(istart+i)] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+1*N[q]+i];
        rc[d+ns_local[rank]*(istart+i)] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+2*N[q]+i];
        rx[d+ns_local[rank]*(istart+i)] = recvbuf[istart*nvar*ns_local[rank]+d*nvar*N[q]+3*N[q]+i];
      }
    }
    istart += N[q];
  }
  /* deallocate receive buffer */
  if (recvbuf)  free(recvbuf);

  /* solve the gathered systems in serial */
  ierr = tridiagLU(ra,rb,rc,rx,NT,ns_local[rank],context,NULL);
  if (ierr) return(ierr);

  /* allocate send buffer and save the data to send */
  if (ns_local[rank] > 0)
    sendbuf = (double*) calloc (ns_local[rank]*NT,sizeof(double));
  else sendbuf = NULL;
  istart = 0;
  for (q = 0; q < nproc; q++) {
    for (i = 0; i < N[q]; i++) {
      for (d = 0; d < ns_local[rank]; d++) {
        sendbuf[istart*ns_local[rank]+d*N[q]+i] = rx[d+ns_local[rank]*(istart+i)];
      }
    }
    istart += N[q];
  }
  dstart = 0;
  for (p = 0; p < nproc; p++) {
    if (ns_local[p] > 0) {

      /* allocate receive buffer */
      recvbuf = (double*) calloc (ns_local[p]*n, sizeof(double));

      /* scatter the solution back */
      for (q = 0; q < nproc; q++) {
        counts[q] = ns_local[p]*N[q];
        displ[q]  = (q == 0 ? 0 : displ[q-1]+counts[q-1]);
      }
      MPI_Scatterv(sendbuf,counts,displ,MPI_DOUBLE,
                   recvbuf,ns_local[p]*n,MPI_DOUBLE,
                   p,*comm);
      /* save the solution on all root processes */
      for (d = 0; d < ns_local[p]; d++) {
        for (i = 0; i < n; i++) {
          x[d+dstart+ns*i] = recvbuf[d*n+i];
        }
      }
      dstart += ns_local[p];
      /* deallocate receive buffer */
      free(recvbuf);
    }
  }
  /* deallocate send buffer */
  if (sendbuf) free(sendbuf);

  /* clean up */
  if (ns_local[rank] > 0) {
    free(ra);
    free(rb);
    free(rc);
    free(rx);
  }
  free(ns_local);
  free(N);
  free(displ);
  free(counts);

  return(0);
#endif
}
