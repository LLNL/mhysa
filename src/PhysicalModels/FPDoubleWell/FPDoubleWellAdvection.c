#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fpdoublewell.h>
#include <hypar.h>

int FPDoubleWellAdvection(double *f,double *u,int dir,void *s,double t)
{
  HyPar         *solver = (HyPar*)        s;
  int           i, v;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double x; _GetCoordinate_(0,index[0]-ghosts,dim,ghosts,solver->x,x);
    for (v = 0; v < nvars; v++) f[nvars*p+v] = drift(x) * u[nvars*p+v];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
