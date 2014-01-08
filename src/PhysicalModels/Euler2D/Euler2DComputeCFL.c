#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

double Euler2DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar   *solver = (HyPar*)   s;
  Euler2D *param  = (Euler2D*) solver->physics;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int index[ndims];
  double *u   = solver->u;

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double rho,vx,vy,e,P,c,dxinv,dyinv,local_cfl[2];
    _Euler2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);

    c     = sqrt(param->gamma*P/rho); /* speed of sound */
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv); /* 1/dy */

    local_cfl[_XDIR_] = (absolute(vx)+c)*dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vy)+c)*dt*dyinv; /* local cfl for this grid point (y) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
