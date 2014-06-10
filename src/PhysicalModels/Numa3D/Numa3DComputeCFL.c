#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

double Numa3DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar  *solver = (HyPar*)  s;
  Numa3D *param  = (Numa3D*) solver->physics;

  int     *dim    = solver->dim_local;
  int     ghosts  = solver->ghosts;
  int     ndims   = solver->ndims;
  double  *u      = solver->u;
  int     index[ndims];

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double drho,uvel,vvel,wvel,dT,rho0,T0,P0,EP,c,zcoord;

    _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->x,zcoord);
    param->StandardAtmosphere(param,zcoord,&EP,&P0,&rho0,&T0);
    _Numa3DGetFlowVars_         ((u+_MODEL_NVARS_*p),drho,uvel,vvel,wvel,dT,rho0);
    _Numa3DComputeSpeedofSound_ (param->gamma,param->R,T0,dT,rho0,drho,EP,c);

    double dxinv, dyinv, dzinv;
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv); /* 1/dy */
    _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->dxinv,dzinv); /* 1/dz */

    double local_cfl[3];
    local_cfl[_XDIR_] = (absolute(uvel)+c)*dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vvel)+c)*dt*dyinv; /* local cfl for this grid point (y) */
    local_cfl[_ZDIR_] = (absolute(wvel)+c)*dt*dzinv; /* local cfl for this grid point (z) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];
    if (local_cfl[_ZDIR_] > max_cfl) max_cfl = local_cfl[_ZDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
