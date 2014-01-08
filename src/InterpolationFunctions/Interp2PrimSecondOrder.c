#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Second order central interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 1

int Interp2PrimSecondOrder(double *fI,double *fC,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p, qL, qR;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      int v; for (v=0; v<nvars; v++) fI[p*nvars+v] = (fC[qR*nvars+v] - fC[qL*nvars+v]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }
  
  return(0);
}
