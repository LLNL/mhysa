#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <secondderivative.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Second order central differencing
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 1

int SecondDerivativeSecondOrderCentral(double *D2f,double *f,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*) s;
  int           i, v;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  if ((!D2f) || (!f)) {
    fprintf(stderr, "Error in SecondDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }
  if (ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in SecondDerivativeSecondOrder(): insufficient number of ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    for (i = 0; i < dim[dir]; i++) {
      int qL, qC, qR;
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC);
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  D2f[qC*nvars+v] = f[qL*nvars+v]-2*f[qC*nvars+v]+f[qR*nvars+v];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }
  
  return(0);
}
