#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/numa2d-cons.h>
#include <hypar.h>

int Numa2DConsFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar       *solver = (HyPar*)   s;
  Numa2DCons  *param  = (Numa2DCons*) solver->physics;
  int         i, p;

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
    _ArrayIndex1DWO_      (ndims,dim,index,offset,ghosts,p);
    _Numa2DConsSetFlux_   ((f+_MODEL_NVARS_*p),dir,param,(u+_MODEL_NVARS_*p));
    _ArrayIncrementIndex_ (ndims,bounds,index,done);
  }
  
  return(0);
}

int Numa2DConsStiffFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar       *solver = (HyPar*)   s;
  Numa2DCons  *param  = (Numa2DCons*) solver->physics;
  int         i, p;

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
    _ArrayIndex1DWO_          (ndims,dim,index,offset,ghosts,p);
    _Numa2DConsSetLinearFlux_ ((f+_MODEL_NVARS_*p),dir,param,(u+_MODEL_NVARS_*p));
    _ArrayIncrementIndex_     (ndims,bounds,index,done);
  }
  
  return(0);
}
