#include <stdio.h>
#include <stdlib.h>
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
  int             ierr    = 0,sys,Nsys,d;

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
  int *indexC       = (int*) calloc (ndims,sizeof(int));
  int *indexI       = (int*) calloc (ndims,sizeof(int));
  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  int *bounds_inter = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  /* calculate total number of tridiagonal systems to solve */
  Nsys = 1; for (d=0; d<ndims; d++) Nsys *= bounds_outer[d]; Nsys *= nvars;

  /* Allocate arrays for tridiagonal system */
  double **A, **B, **C, **R;
  A = (double**) calloc (Nsys,sizeof(double*));
  B = (double**) calloc (Nsys,sizeof(double*));
  C = (double**) calloc (Nsys,sizeof(double*));
  R = (double**) calloc (Nsys,sizeof(double*));
  for (sys = 0; sys < Nsys; sys++) {
    A[sys] = (double*) calloc (dim[dir]+1, sizeof(double));
    B[sys] = (double*) calloc (dim[dir]+1, sizeof(double));
    C[sys] = (double*) calloc (dim[dir]+1, sizeof(double));
    R[sys] = (double*) calloc (dim[dir]+1, sizeof(double));
  }

  int done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  sys = 0;
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,indexC,ndims); CHECKERR(ierr);
    ierr = ArrayCopy1D_int(index_outer,indexI,ndims); CHECKERR(ierr);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; qm3 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-2; qm2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-1; qm1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]  ; qp1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]+1; qp2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      } else {
        indexC[dir] = indexI[dir]+2; qm3 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]+1; qm2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]  ; qm1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-1; qp1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-2; qp2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
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
          A[sys*nvars+v][indexI[dir]] = 0.0;
          B[sys*nvars+v][indexI[dir]] = 1.0;
          C[sys*nvars+v][indexI[dir]] = 0.0;
        } else {
          A[sys*nvars+v][indexI[dir]] = (2*one_third)*w1 + (one_third)*w2;
          B[sys*nvars+v][indexI[dir]] = (one_third)*w1 + (2*one_third)*(w2+w3);
          C[sys*nvars+v][indexI[dir]] = (one_third)*w3;
        }
        R[sys*nvars+v][indexI[dir]] = w1*f1 + w2*f2 + w3*f3;
      }
    }
    sys++;
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

#ifdef serial

  /* Solve the tridiagonal system */
  ierr = tridiagLU(A,B,C,R,dim[dir]+1,Nsys,NULL,NULL);

#else

  /* Set the MPI context for the tridiagonal system solver */
  MPI_Comm world; MPI_Comm_dup(MPI_COMM_WORLD,&world);

  MPIContext mpicntxt;
  mpicntxt.rank  = mpi->ip[dir];     /* rank along this dimension  */
  mpicntxt.nproc = mpi->iproc[dir];  /* nproc along this dimension */
  mpicntxt.comm  = &world;

  int *ip = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(mpi->ip,ip,ndims); CHECKERR(ierr);
  mpicntxt.proc  = (int*) calloc (mpicntxt.nproc,sizeof(int));
  for (d=0; d<mpicntxt.nproc; d++) {
    ip[dir] = d;
    int rank = MPIRank1D(ndims,mpi->iproc,ip);
    mpicntxt.proc[d] = rank;
  }

  /* Solve the tridiagonal system */
  /* all processes except the last will solve without the last interface to avoid overlap */
  if (mpi->ip[dir] != mpi->iproc[dir]-1)  ierr = tridiagLU(A,B,C,R,dim[dir]  ,Nsys,NULL,&mpicntxt);
  else                                    ierr = tridiagLU(A,B,C,R,dim[dir]+1,Nsys,NULL,&mpicntxt);
  /* Now get the solution to the last interface from the next proc */
  ierr = ArrayCopy1D_int(mpi->ip,ip,ndims); CHECKERR(ierr);
  ip[dir]++; int source = MPIRank1D(ndims,mpi->iproc,ip); ip[dir]--;
  ip[dir]--; int dest   = MPIRank1D(ndims,mpi->iproc,ip); ip[dir]++;
  double *sendbuf,*recvbuf;
  sendbuf = (double*) calloc (Nsys,sizeof(double));
  recvbuf = (double*) calloc (Nsys,sizeof(double));
  MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
  if (mpi->ip[dir]) for (d=0; d<Nsys; d++) sendbuf[d] = R[d][0];
  if (mpi->ip[dir] != mpi->iproc[dir]-1) MPI_Irecv(recvbuf,Nsys,MPI_DOUBLE,source,214,MPI_COMM_WORLD,&req[0]);
  if (mpi->ip[dir])                      MPI_Isend(sendbuf,Nsys,MPI_DOUBLE,dest  ,214,MPI_COMM_WORLD,&req[1]);
  MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  if (mpi->ip[dir] != mpi->iproc[dir]-1) for (d=0; d<Nsys; d++) R[d][dim[dir]] = recvbuf[d];
  free(sendbuf);
  free(recvbuf);

  /* deallocate allocations made for MPI context */
  free(ip);
  free(mpicntxt.proc);
  MPI_Comm_free(&world);

#endif

  if (ierr == -1) {
    /* check if singular */
    fprintf(stderr,"Error in Interp1PrimFifthOrderCRWENO(): Singular system!\n");
    return(ierr);
  } else if (ierr) {
    /* if any other error */
    return(ierr);
  } else {
    /* save the solution to fI */
    int done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
    sys = 0;
    while (!done) {
      ierr = ArrayCopy1D_int(index_outer,indexI,ndims); CHECKERR(ierr);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p = ArrayIndex1D(ndims,bounds_inter,indexI,NULL,0);
        int v; for (v=0; v<nvars; v++) fI[nvars*p+v] = R[sys*nvars+v][indexI[dir]];
      }
      done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
      sys++;
    }
  }

  for (sys = 0; sys < Nsys; sys++) {
    free(A[sys]);
    free(B[sys]);
    free(C[sys]);
    free(R[sys]);
  }
  free(A);
  free(B);
  free(C);
  free(R);

  free(indexC);
  free(indexI);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);
  
  return(0);
}
