#include <stdlib.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;
  int       i;

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
    double rho, v, e, P;
    _Euler1DGetFlowVar_((u+_MODEL_NVARS_*p),rho,v,e,P,param);
    _Euler1DSetFlux_((f+_MODEL_NVARS_*p),rho,v,e,P);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

int Euler1DStiffFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;
  int       i;

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
    double rho, v, e, P;
    _Euler1DGetFlowVar_((u+_MODEL_NVARS_*p),rho,v,e,P,param);
    _Euler1DSetStiffFlux_((f+_MODEL_NVARS_*p),rho,v,e,P);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
