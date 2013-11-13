#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

int Euler2DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar   *solver = (HyPar*)   s;
  Euler2D *param  = (Euler2D*) solver->physics;
  int     i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;

  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, vx, vy, e, P;
    _Euler2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);
    _Euler2DSetFlux_((f+_MODEL_NVARS_*p),rho,vx,vy,e,P,param,dir);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
