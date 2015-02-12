#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DPreStep(double *u,void *s,void *m,double waqt)
{
  HyPar             *solver = (HyPar*)   s;
  Euler1D           *param  = (Euler1D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);
  /* copy the solution to act as a reference for linearization */
  _ArrayCopy1D_(u,param->solution,(solver->npoints_local_wghosts*_MODEL_NVARS_));

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, v, e, P;
    _Euler1DGetFlowVar_((u+_MODEL_NVARS_*p),rho,v,e,P,param);
    _Euler1DSetStiffJac_((param->fast_jac+JacSize*p),rho,v,e,P,param->gamma);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
