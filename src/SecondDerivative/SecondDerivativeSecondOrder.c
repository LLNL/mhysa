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

int SecondDerivativeSecondOrder(double *D2f,double *f,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  int           ierr    = 0, i, v;

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
  int *indexC       = (int*) calloc (ndims,sizeof(int));
  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;

  int done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,indexC,ndims); CHECKERR(ierr);
    for (i = 0; i < dim[dir]; i++) {
      int p, qL, qC, qR;
      indexC[dir] = i  ; p  = ArrayIndex1D(ndims,dim,indexC,NULL,0     );
      indexC[dir] = i-1; qL = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      indexC[dir] = i  ; qC = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      indexC[dir] = i+1; qR = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      for (v=0; v<nvars; v++)  D2f[p*nvars+v] = f[qL*nvars+v]-2*f[qC*nvars+v]+f[qR*nvars+v];
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  free(indexC);
  free(index_outer);
  free(bounds_outer);
  
  return(0);
}
