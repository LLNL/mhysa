#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

int Numa3DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar  *solver = (HyPar*)   s;
  Numa3D *param  = (Numa3D*) solver->physics;
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
    double drho,uvel,vvel,wvel,dT,dP,rho0,T0,P0,EP,zcoord;

    _GetCoordinate_(_ZDIR_,index[_ZDIR_]-ghosts,dim,ghosts,solver->x,zcoord);
    param->StandardAtmosphere(param,zcoord,&EP,&P0,&rho0,&T0);

    _Numa3DGetFlowVars_     ((u+_MODEL_NVARS_*p),drho,uvel,vvel,wvel,dT,rho0);
    _Numa3DComputePressure_ (param,T0,dT,P0,dP);
    _Numa3DSetFlux_         ((f+_MODEL_NVARS_*p),dir,drho,uvel,vvel,wvel,dT,dP,rho0,T0);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  
  return(0);
}

int Numa3DStiffFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar  *solver = (HyPar*)   s;
  Numa3D *param  = (Numa3D*) solver->physics;
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
    double drho,uvel,vvel,wvel,dT,dP,rho0,T0,P0,EP,zcoord;

    _GetCoordinate_(_ZDIR_,index[_ZDIR_]-ghosts,dim,ghosts,solver->x,zcoord);
    param->StandardAtmosphere(param,zcoord,&EP,&P0,&rho0,&T0);

    _Numa3DGetFlowVars_               ((u+_MODEL_NVARS_*p),drho,uvel,vvel,wvel,dT,rho0);
    _Numa3DComputeLinearizedPressure_ (param,T0,dT,P0,dP);
    _Numa3DSetLinearFlux_             ((f+_MODEL_NVARS_*p),dir,drho,uvel,vvel,wvel,dT,dP,rho0,T0);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  
  return(0);
}
