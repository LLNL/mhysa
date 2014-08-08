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
  First order upwind characteristic-based interpolation (uniform grid )
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 1

int Interp1PrimFirstOrderUpwindChar(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*) s;
  int           i, k, v;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double R[nvars*nvars], L[nvars*nvars], uavg[nvars], fchar[nvars];

#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      int pL, pR; /* 1D index of the left and right cells */
      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,pL);
      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,pR);
      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&u[nvars*pL],&u[nvars*pR],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double fcL = 0, fcR = 0;
        for (k = 0; k < nvars; k++) {
          fcL += L[v*nvars+k] * fC[pL*nvars+k];
          fcR += L[v*nvars+k] * fC[pR*nvars+k];
        }

        /* first order upwind approximation of the characteristic flux */
        fchar[v] = (upw > 0 ? fcL : fcR);

      }

      /* calculate the interface u from the characteristic u */
      IERR MatVecMult(nvars,(fI+nvars*p),R,fchar); CHECKERR(ierr);

    }
  }

  return(0);
}
