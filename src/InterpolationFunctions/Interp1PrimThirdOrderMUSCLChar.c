#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_omp
#include <omp.h>
#endif

/* 
  Third order MUSCL characteristic-based interpolation (uniform grid )
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int Interp1PrimThirdOrderMUSCLChar(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m,int uflag)
{
  HyPar             *solver = (HyPar*) s;
  MUSCLParameters   *muscl   = (MUSCLParameters*) solver->interp;
  int               i, k, v;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  double one_third = 1.0/3.0;
  double one_sixth = 1.0/6.0;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double R[nvars*nvars], L[nvars*nvars], uavg[nvars], fchar[nvars];

  if (upw > 0) {
#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
    for (i=0; i<N_outer; i++) {
      _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
      _ArrayCopy1D_(index_outer,indexC,ndims);
      _ArrayCopy1D_(index_outer,indexI,ndims);

      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

        /* 1D indices of the stencil grid points */
        int qm1,qm2,qp1;
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);

        int p; /* 1D index of the interface */
        _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

        /* find averaged state at this interface */
        IERR solver->AveragingFunction(uavg,&u[nvars*qm1],&u[nvars*qp1],solver->physics);
        CHECKERR(ierr);

        /* Get the left and right eigenvectors */
        IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
        IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

        /* For each characteristic field */
        for (v = 0; v < nvars; v++) {
          /* calculate the characteristic flux components along this characteristic */
          double m2, m1, p1;
          m2 = m1 = p1 = 0;
          for (k = 0; k < nvars; k++) {
            m2 += L[v*nvars+k] * fC[qm2*nvars+k];
            m1 += L[v*nvars+k] * fC[qm1*nvars+k];
            p1 += L[v*nvars+k] * fC[qp1*nvars+k];
          }
          double fdiff = p1 - m1;
          double bdiff = m1 - m2;
          double limit =  (3*fdiff*bdiff + muscl->eps) 
                        / (2*(fdiff-bdiff)*(fdiff-bdiff) + 3*fdiff*bdiff + muscl->eps);
          fchar[v] = m1 +  limit * (one_third*fdiff + one_sixth*bdiff);
        }

        /* calculate the interface u from the characteristic u */
        IERR MatVecMult(nvars,(fI+nvars*p),R,fchar); CHECKERR(ierr);
      }
    }
  } else {
#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
    for (i=0; i<N_outer; i++) {
      _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
      _ArrayCopy1D_(index_outer,indexC,ndims);
      _ArrayCopy1D_(index_outer,indexI,ndims);

      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

        /* 1D indices of the stencil grid points */
        int qm1,qp1,qp2;
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);

        int p; /* 1D index of the interface */
        _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

        /* find averaged state at this interface */
        IERR solver->AveragingFunction(uavg,&u[nvars*qm1],&u[nvars*qp1],solver->physics);
        CHECKERR(ierr);

        /* Get the left and right eigenvectors */
        IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
        IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

        /* For each characteristic field */
        for (v = 0; v < nvars; v++) {
          /* calculate the characteristic flux components along this characteristic */
          double m1, p1, p2;
          m1 = p1 = p2 = 0;
          for (k = 0; k < nvars; k++) {
            m1 += L[v*nvars+k] * fC[qm1*nvars+k];
            p1 += L[v*nvars+k] * fC[qp1*nvars+k];
            p2 += L[v*nvars+k] * fC[qp2*nvars+k];
          }
          double fdiff = p2 - p1;
          double bdiff = p1 - m1;
          double limit =  (3*fdiff*bdiff + muscl->eps) 
                        / (2*(fdiff-bdiff)*(fdiff-bdiff) + 3*fdiff*bdiff + muscl->eps);
          fchar[v] = p1 -  limit * (one_third*fdiff + one_sixth*bdiff);
        }

        /* calculate the interface u from the characteristic u */
        IERR MatVecMult(nvars,(fI+nvars*p),R,fchar); CHECKERR(ierr);
      }
    }
  }

  return(0);
}
