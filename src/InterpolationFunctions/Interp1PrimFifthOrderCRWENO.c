#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <tridiagLU.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_omp
#include <omp.h>
#endif

/* 
  Fifth order CRWENO interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int Interp1PrimFifthOrderCRWENO(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  TridiagLU       *lu     = (TridiagLU*)      solver->lusolver;
  int             sys,Nsys,d;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  static const double one_third          = 1.0/3.0;
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  double *ww1, *ww2, *ww3;
  ww1 = weno->w1 + (upw < 0 ? 2*weno->size : 0) + (u==fC ? weno->size : 0) + weno->offset[dir];
  ww2 = weno->w2 + (upw < 0 ? 2*weno->size : 0) + (u==fC ? weno->size : 0) + weno->offset[dir];
  ww3 = weno->w3 + (upw < 0 ? 2*weno->size : 0) + (u==fC ? weno->size : 0) + weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  /* calculate total number of tridiagonal systems to solve */
  _ArrayProduct1D_(bounds_outer,ndims,Nsys); Nsys *= nvars;

  /* Allocate arrays for tridiagonal system */
  double *A = weno->A;
  double *B = weno->B;
  double *C = weno->C;
  double *R = weno->R;

#pragma omp parallel for schedule(auto) default(shared) private(sys,d,index_outer,indexC,indexI)
  for (sys=0; sys < N_outer; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
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

      /* Defining stencil points */
      double *fm3, *fm2, *fm1, *fp1, *fp2;
      fm3 = fC+qm3*nvars;
      fm2 = fC+qm2*nvars;
      fm1 = fC+qm1*nvars;
      fp1 = fC+qp1*nvars;
      fp2 = fC+qp2*nvars;

      /* Candidate stencils and their optimal weights*/
      double f1[nvars], f2[nvars], f3[nvars];
      if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
          || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
        /* Use WENO5 at the physical boundaries */
        _ArrayAXBYCZ_(f1,(2*one_sixth),fm3,(-7*one_sixth) ,fm2,(11*one_sixth) ,fm1,nvars);
        _ArrayAXBYCZ_(f2,(-one_sixth) ,fm2,(5*one_sixth)  ,fm1,(2*one_sixth)  ,fp1,nvars);
        _ArrayAXBYCZ_(f3,(2*one_sixth),fm1,(5*one_sixth)  ,fp1,(-one_sixth)   ,fp2,nvars);
      } else {
        /* CRWENO5 at the interior points */
        _ArrayAXBY_(f1,(one_sixth)  ,fm2,(5*one_sixth),fm1,nvars);
        _ArrayAXBY_(f2,(5*one_sixth),fm1,(one_sixth)  ,fp1,nvars);
        _ArrayAXBY_(f3,(one_sixth)  ,fm1,(5*one_sixth),fp1,nvars);
      }

      /* calculate WENO weights */
      double *w1, *w2, *w3;
      w1 = (ww1+p*nvars);
      w2 = (ww2+p*nvars);
      w3 = (ww3+p*nvars);

      if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
          || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
        _ArraySetValue_     ((A+Nsys*indexI[dir]+sys*nvars),nvars,0.0)
        _ArraySetValue_     ((B+Nsys*indexI[dir]+sys*nvars),nvars,1.0)
        _ArraySetValue_     ((C+Nsys*indexI[dir]+sys*nvars),nvars,0.0)
      } else {
        if (upw > 0) {
          _ArrayAXBY_       ((A+Nsys*indexI[dir]+sys*nvars),(2*one_third) ,w1,(one_third)  ,w2,nvars);
          _ArrayAXBYCZ_     ((B+Nsys*indexI[dir]+sys*nvars),(one_third)   ,w1,(2*one_third),w2,(2*one_third),w3,nvars);
          _ArrayScaleCopy1D_(w3,(one_third),(C+Nsys*indexI[dir]+sys*nvars),nvars);
        } else {
          _ArrayAXBY_       ((C+Nsys*indexI[dir]+sys*nvars),(2*one_third) ,w1,(one_third)  ,w2,nvars);
          _ArrayAXBYCZ_     ((B+Nsys*indexI[dir]+sys*nvars),(one_third)   ,w1,(2*one_third),w2,(2*one_third),w3,nvars);
          _ArrayScaleCopy1D_(w3,(one_third),(A+Nsys*indexI[dir]+sys*nvars),nvars);
        }
      }
      _ArrayMultiply3Add1D_ ((R+Nsys*indexI[dir]+sys*nvars),w1,f1,w2,f2,w3,f3,nvars);
    }
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
#pragma omp parallel for schedule(auto) default(shared) private(sys,d,index_outer,indexC,indexI)
  for (sys=0; sys < N_outer; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      _ArrayCopy1D_((R+sys*nvars+Nsys*indexI[dir]),(fI+nvars*p),nvars);
    }
  }

  return(0);
}
