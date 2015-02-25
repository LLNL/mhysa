#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_omp
#include <omp.h>
#endif

/* 
  Second order central interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 1

int Interp1PrimSecondOrderCentral(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m,int uflag)
{
  HyPar         *solver = (HyPar*)        s;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexL[ndims], indexR[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  int i;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexL,indexR,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexL,ndims);
    _ArrayCopy1D_(index_outer,indexR,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      indexL[dir] = indexI[dir]-1;
      indexR[dir] = indexI[dir];
      int p;  _ArrayIndex1D_(ndims,bounds_inter,indexI,0     ,p);
      int qL; _ArrayIndex1D_(ndims,dim         ,indexL,ghosts,qL);
      int qR; _ArrayIndex1D_(ndims,dim         ,indexR,ghosts,qR);
      int v; for (v=0; v<nvars; v++)  fI[p*nvars+v] = 0.5*(fC[qL*nvars+v]+fC[qR*nvars+v]);
    }
  }

  return(0);
}
