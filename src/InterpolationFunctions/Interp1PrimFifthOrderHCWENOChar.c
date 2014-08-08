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
  Fifth order HCWENO characteristic-based interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int Interp1PrimFifthOrderHCWENOChar(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  TridiagLU       *lu     = (TridiagLU*)      solver->lusolver;
  int             sys,Nsys,d,v,k;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  static const double one_half           = 1.0/2.0;
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

  /* calculate total number of block tridiagonal systems to solve */
  _ArrayProduct1D_(bounds_outer,ndims,Nsys);

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double R[nvars*nvars], L[nvars*nvars], uavg[nvars], fchar[nvars];

  /* Allocate arrays for tridiagonal system */
  double *A = weno->A;
  double *B = weno->B;
  double *C = weno->C;
  double *F = weno->R;

#pragma omp parallel for schedule(auto) default(shared) private(sys,d,v,k,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (sys=0; sys<Nsys; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
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

      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&u[nvars*qm1],&u[nvars*qp1],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      for (v=0; v<nvars; v++)  {

        /* calculate the characteristic flux components along this characteristic */
        double fm3, fm2, fm1, fp1, fp2;
        fm3 = fm2 = fm1 = fp1 = fp2 = 0;
        for (k = 0; k < nvars; k++) {
          fm3 += L[v*nvars+k] * fC[qm3*nvars+k];
          fm2 += L[v*nvars+k] * fC[qm2*nvars+k];
          fm1 += L[v*nvars+k] * fC[qm1*nvars+k];
          fp1 += L[v*nvars+k] * fC[qp1*nvars+k];
          fp2 += L[v*nvars+k] * fC[qp2*nvars+k];
        }

        /* Candidate stencils and their optimal weights*/
        double f1, f2, f3;
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          f1 = (2*one_sixth)*fm3 - (7.0*one_sixth)*fm2 + (11.0*one_sixth)*fm1;
          f2 = (-one_sixth)*fm2 + (5.0*one_sixth)*fm1 + (2*one_sixth)*fp1;
          f3 = (2*one_sixth)*fm1 + (5*one_sixth)*fp1 - (one_sixth)*fp2;
        } else {
          /* HCWENO5 at the interior points */
          f1 = (one_sixth) * (fm2 + 5*fm1);
          f2 = (one_sixth) * (5*fm1 + fp1);
          f3 = (one_sixth) * (fm1 + 5*fp1);
        }

        /* calculate WENO weights */
        double w1,w2,w3;
        w1 = *(ww1+p*nvars+v);
        w2 = *(ww2+p*nvars+v);
        w3 = *(ww3+p*nvars+v);

        /* calculate the hybridization parameter */
        double sigma;
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Standard WENO at physical boundaries */
          sigma = 0.0;
        } else {
          double cuckoo, df_jm12, df_jp12, df_jp32, r_j, r_jp1, r_int;
          cuckoo = (0.9*weno->rc / (1.0-0.9*weno->rc)) * weno->xi * weno->xi;
          df_jm12 = fm1 - fm2;
          df_jp12 = fp1 - fm1;
          df_jp32 = fp2 - fp1;
          r_j   = (absolute(2*df_jp12*df_jm12)+cuckoo)/(df_jp12*df_jp12+df_jm12*df_jm12+cuckoo);
          r_jp1 = (absolute(2*df_jp32*df_jp12)+cuckoo)/(df_jp32*df_jp32+df_jp12*df_jp12+cuckoo);
          r_int = min(r_j, r_jp1);
          sigma = min((r_int/weno->rc), 1.0); 
        }

        if (upw > 0) {
          for (k=0; k<nvars; k++) {
            A[(Nsys*indexI[dir]+sys)*nvars*nvars+v*nvars+k] = (one_half*sigma)  * L[v*nvars+k];
            B[(Nsys*indexI[dir]+sys)*nvars*nvars+v*nvars+k] = (1.0)             * L[v*nvars+k];
            C[(Nsys*indexI[dir]+sys)*nvars*nvars+v*nvars+k] = (one_sixth*sigma) * L[v*nvars+k];
          }
        } else {
          for (k=0; k<nvars; k++) {
            C[(Nsys*indexI[dir]+sys)*nvars*nvars+v*nvars+k] = (one_half*sigma)  * L[v*nvars+k];
            B[(Nsys*indexI[dir]+sys)*nvars*nvars+v*nvars+k] = (1.0)             * L[v*nvars+k];
            A[(Nsys*indexI[dir]+sys)*nvars*nvars+v*nvars+k] = (one_sixth*sigma) * L[v*nvars+k];
          }
        }
        double fWENO, fCompact;
        fWENO    = w1*f1 + w2*f2 + w3*f3;
        fCompact = one_sixth * (one_third*fm2 + 19.0*one_third*fm1 + 10.0*one_third*fp1);
        F[(Nsys*indexI[dir]+sys)*nvars+v] = sigma*fCompact + (1.0-sigma)*fWENO;
      }
    }
  }

#ifdef serial

  /* Solve the tridiagonal system */
  IERR blocktridiagLU(A,B,C,F,dim[dir]+1,Nsys,nvars,lu,NULL); CHECKERR(ierr);

#else

  /* Solve the tridiagonal system */
  /* all processes except the last will solve without the last interface to avoid overlap */
  if (mpi->ip[dir] != mpi->iproc[dir]-1)  { 
    IERR blocktridiagLU(A,B,C,F,dim[dir]  ,Nsys,nvars,lu,&mpi->comm[dir]); CHECKERR(ierr); 
  } else { 
    IERR blocktridiagLU(A,B,C,F,dim[dir]+1,Nsys,nvars,lu,&mpi->comm[dir]); CHECKERR(ierr);
  }

  /* Now get the solution to the last interface from the next proc */
  double *sendbuf = weno->sendbuf;
  double *recvbuf = weno->recvbuf;
  MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
  if (mpi->ip[dir]) for (d=0; d<Nsys*nvars; d++) sendbuf[d] = F[d];
  if (mpi->ip[dir] != mpi->iproc[dir]-1) MPI_Irecv(recvbuf,Nsys*nvars,MPI_DOUBLE,mpi->ip[dir]+1,214,mpi->comm[dir],&req[0]);
  if (mpi->ip[dir])                      MPI_Isend(sendbuf,Nsys*nvars,MPI_DOUBLE,mpi->ip[dir]-1,214,mpi->comm[dir],&req[1]);
  MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  if (mpi->ip[dir] != mpi->iproc[dir]-1) for (d=0; d<Nsys*nvars; d++) F[d+Nsys*nvars*dim[dir]] = recvbuf[d];

#endif

  /* save the solution to fI */
#pragma omp parallel for schedule(auto) default(shared) private(sys,d,v,k,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (sys=0; sys<Nsys; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      int v; for (v=0; v<nvars; v++) fI[nvars*p+v] = F[sys*nvars+v+Nsys*nvars*indexI[dir]];
    }
  }

  return(0);
}
