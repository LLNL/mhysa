#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

int FPPowerSystem1BusAdvection(double *f,double *u,int dir,void *s,double t)
{
  HyPar *solver = (HyPar*) s;
  int   *dim    = solver->dim_local;
  int   ghosts  = solver->ghosts;
  int   ndims   = solver->ndims;
  int   nvars   = solver->nvars;

  /* calculate total size of arrays */
  int bounds[ndims]; _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  int size;          _ArrayProduct1D_(bounds,ndims,size); size *= nvars;

  /* f = u */
  _ArrayCopy1D_(u,f,size);

  return(0);
}
