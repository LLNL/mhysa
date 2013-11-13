#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

double FPPowerSystemDriftFunction(int,void*,double,double,double);

int FPPowerSystemAdvection(double *f,double *u,int dir,void *s,double t)
{
  HyPar         *solver = (HyPar*)        s;
//  FPPowerSystem *params = (FPPowerSystem*)solver->physics;
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
//    double x     = solver->GetCoordinate(0,index[0],dim,ghosts,solver->x);
//    double y     = solver->GetCoordinate(1,index[1],dim,ghosts,solver->x);
//    double drift = FPPowerSystemDriftFunction(dir,params,x,y,t);
//    for (v = 0; v < nvars; v++) f[nvars*p+v] = drift * u[nvars*p+v];
    for (v = 0; v < nvars; v++) f[nvars*p+v] = u[nvars*p+v];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
