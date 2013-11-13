#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <tridiagLU.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Fifth order CRWENO interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int Interp1PrimFifthOrderCRWENO(double *fI,double *fC,double *u,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  TridiagLU       *lu     = (TridiagLU*)      solver->lusolver;
  int             sys,Nsys,d,done;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  double one_third          = 1.0/3.0;
  double one_sixth          = 1.0/6.0;
  double thirteen_by_twelve = 13.0/12.0;
  double one_fourth         = 1.0/4.0;

  if ((!fI) || (!fC)) {
    fprintf(stderr, "Error in Interp1PrimFifthOrderWENO(): input arrays not allocated.\n");
    return(1);
  }
  if (ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in Interp1PrimFifthOrderWENO(): insufficient number of ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  /* calculate total number of tridiagonal systems to solve */
  Nsys = 1; for (d=0; d<ndims; d++) Nsys *= bounds_outer[d]; Nsys *= nvars;

  /* Allocate arrays for tridiagonal system */
  double *A = weno->A;
  double *B = weno->B;
  double *C = weno->C;
  double *R = weno->R;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  sys = 0;
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }
      int v; 
      for (v=0; v<nvars; v++)  {
        /* Defining stencil points */
        double m3, m2, m1, p1, p2;
        m3 = fC[qm3*nvars+v];
        m2 = fC[qm2*nvars+v];
        m1 = fC[qm1*nvars+v];
        p1 = fC[qp1*nvars+v];
        p2 = fC[qp2*nvars+v];

        /* Candidate stencils and their optimal weights*/
        double f1, f2, f3, c1, c2, c3;
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          f1 = (2*one_sixth)*m3 - (7.0*one_sixth)*m2 + (11.0*one_sixth)*m1; c1 = 0.1;
          f2 = (-one_sixth)*m2 + (5.0*one_sixth)*m1 + (2*one_sixth)*p1;     c2 = 0.6;
          f3 = (2*one_sixth)*m1 + (5*one_sixth)*p1 - (one_sixth)*p2;        c3 = 0.3;
        } else {
          /* CRWENO5 at the interior points */
          f1 = (one_sixth) * (m2 + 5*m1);                                   c1 = 0.2;
          f2 = (one_sixth) * (5*m1 + p1);                                   c2 = 0.5;
          f3 = (one_sixth) * (m1 + 5*p1);                                   c3 = 0.3;
        }

        /* calculate WENO weights */
        double w1,w2,w3;

        if (weno->no_limiting) {
          /* fifth order polynomial interpolation */
          w1 = c1;
          w2 = c2;
          w3 = c3;
        } else {
          /* Smoothness indicators */
          double b1, b2, b3;
          b1 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
               + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
          b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
               + one_fourth*(m2-p1)*(m2-p1);
          b3 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
               + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
  
          /* This parameter is needed for Borges' or Yamaleev-Carpenter
             implementation of non-linear weights  */
          double tau;
          if (weno->borges) {
            tau = absolute(b3 - b1);
          } else if (weno->yc) {
            tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
          } else {
            tau = 0;
          }
  
          /* Defining the non-linear weights */
          double a1, a2, a3;
          if (weno->borges || weno->yc) {
            a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));
            a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));
            a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));
          } else {
            a1 = c1 / raiseto(b1+weno->eps,weno->p);
            a2 = c2 / raiseto(b2+weno->eps,weno->p);
            a3 = c3 / raiseto(b3+weno->eps,weno->p);
          }

          /* Convexity */
          double a_sum_inv;
          a_sum_inv = 1.0 / (a1 + a2 + a3);
          w1 = a1 * a_sum_inv;
          w2 = a2 * a_sum_inv;
          w3 = a3 * a_sum_inv;
  
          /* Mapping the weights, if needed */
          if (weno->mapped) {
            a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
            a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
            a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
            a_sum_inv = 1.0 / (a1 + a2 + a3);
            w1 = a1 * a_sum_inv;
            w2 = a2 * a_sum_inv;
            w3 = a3 * a_sum_inv;
          }
  
        }

        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          A[sys*nvars+v+Nsys*indexI[dir]] = 0.0;
          B[sys*nvars+v+Nsys*indexI[dir]] = 1.0;
          C[sys*nvars+v+Nsys*indexI[dir]] = 0.0;
        } else {
          if (upw > 0) {
            A[sys*nvars+v+Nsys*indexI[dir]] = (2*one_third)*w1 + (one_third)*w2;
            B[sys*nvars+v+Nsys*indexI[dir]] = (one_third)*w1 + (2*one_third)*(w2+w3);
            C[sys*nvars+v+Nsys*indexI[dir]] = (one_third)*w3;
          } else {
            C[sys*nvars+v+Nsys*indexI[dir]] = (2*one_third)*w1 + (one_third)*w2;
            B[sys*nvars+v+Nsys*indexI[dir]] = (one_third)*w1 + (2*one_third)*(w2+w3);
            A[sys*nvars+v+Nsys*indexI[dir]] = (one_third)*w3;
          }
        }
        R[sys*nvars+v+Nsys*indexI[dir]] = w1*f1 + w2*f2 + w3*f3;
      }
    }
    sys++;
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

#ifdef serial

  /* Solve the tridiagonal system */
  IERR tridiagLU(A,B,C,R,dim[dir]+1,Nsys,lu,NULL); CHECKERR(ierr);

#else

  /* Solve the tridiagonal system */
  /* all processes except the last will solve without the last interface to avoid overlap */
  if (mpi->ip[dir] != mpi->iproc[dir]-1)  { IERR tridiagLU(A,B,C,R,dim[dir]  ,Nsys,lu,&mpi->comm[dir]); CHECKERR(ierr); }
  else                                    { IERR tridiagLU(A,B,C,R,dim[dir]+1,Nsys,lu,&mpi->comm[dir]); CHECKERR(ierr); }

  /* Now get the solution to the last interface from the next proc */
  double *sendbuf = weno->sendbuf;
  double *recvbuf = weno->recvbuf;
  MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
  if (mpi->ip[dir]) for (d=0; d<Nsys; d++) sendbuf[d] = R[d];
  if (mpi->ip[dir] != mpi->iproc[dir]-1) MPI_Irecv(recvbuf,Nsys,MPI_DOUBLE,mpi->ip[dir]+1,214,mpi->comm[dir],&req[0]);
  if (mpi->ip[dir])                      MPI_Isend(sendbuf,Nsys,MPI_DOUBLE,mpi->ip[dir]-1,214,mpi->comm[dir],&req[1]);
  MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  if (mpi->ip[dir] != mpi->iproc[dir]-1) for (d=0; d<Nsys; d++) R[d+Nsys*dim[dir]] = recvbuf[d];

#endif

  /* save the solution to fI */
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  sys = 0;
  while (!done) {
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      int v; for (v=0; v<nvars; v++) fI[nvars*p+v] = R[sys*nvars+v+Nsys*indexI[dir]];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
    sys++;
  }

  return(0);
}
