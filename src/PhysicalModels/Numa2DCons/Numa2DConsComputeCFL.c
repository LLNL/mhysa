#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d-cons.h>
#include <hypar.h>

double Numa2DConsComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar       *solver = (HyPar*)  s;
  Numa2DCons  *param  = (Numa2DCons*) solver->physics;

  int     *dim    = solver->dim_local;
  int     ghosts  = solver->ghosts;
  int     ndims   = solver->ndims;
  double  *u      = solver->u;
  int     index[ndims];

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double rho,uvel,vvel,theta,c;

    _Numa2DConsGetFlowVars_         ((u+_MODEL_NVARS_*p),rho,uvel,vvel,theta);
    _Numa2DConsComputeSpeedofSound_ (param,rho,theta,c);

    double dxinv, dyinv;
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv); /* 1/dy */

    double local_cfl[2];
    local_cfl[_XDIR_] = (absolute(uvel)+c)*dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vvel)+c)*dt*dyinv; /* local cfl for this grid point (y) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
