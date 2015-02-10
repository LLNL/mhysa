#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

double FPPowerSystem1BusDissipationFunction(int,int,void*,double);

int FPPowerSystem1BusDiffusionLaplacian(double *f,double *u,int dir,void *s,double t)
{
  HyPar             *solver = (HyPar*)             s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) solver->physics;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;

  /* calculate total size of arrays */
  int bounds[ndims]; _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  int size;          _ArrayProduct1D_(bounds,ndims,size); size *= nvars;

  /* calculate dissipation coefficient  -- constant in x and y */
  double dissipation = FPPowerSystem1BusDissipationFunction(dir,dir,params,t);

  /* f = dissipation * u */
  _ArrayScaleCopy1D_(u,dissipation,f,size);

  return(0);
}

int FPPowerSystem1BusDiffusionGeneral(double *f,double *u,int dir1,int dir2,void *s,double t)
{
  HyPar             *solver = (HyPar*)             s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) solver->physics;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;

  /* calculate total size of arrays */
  int bounds[ndims]; _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  int size;          _ArrayProduct1D_(bounds,ndims,size); size *= nvars;

  /* calculate dissipation coefficient  -- constant in x and y */
  double dissipation = FPPowerSystem1BusDissipationFunction(dir1,dir2,params,t);

  /* f = dissipation * u */
  _ArrayScaleCopy1D_(u,dissipation,f,size);

  return(0);
}
