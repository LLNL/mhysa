#include <stdlib.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar             *solver = (HyPar*)   s;
  Euler1D           *param  = (Euler1D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  nvars   = _MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, v, e, P;
    _Euler1DGetFlowVar_((u+nvars*p),rho,v,e,P,param);
    _Euler1DSetFlux_((f+nvars*p),rho,v,e,P);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

int Euler1DStiffFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar             *solver = (HyPar*)   s;
  Euler1D           *param  = (Euler1D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  nvars   = _MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
//    double rho, v, e, P;
//    _Euler1DGetFlowVar_((u+nvars*p),rho,v,e,P,param);
//    _Euler1DSetStiffFlux_((f+nvars*p),rho,v,e,P,param->gamma);
    _Euler1DSetLinearizedStiffFlux_((f+nvars*p),(u+nvars*p),(param->fast_jac+nvars*nvars*p));
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
