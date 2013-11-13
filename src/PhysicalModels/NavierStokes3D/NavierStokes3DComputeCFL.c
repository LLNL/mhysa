#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

int NavierStokes3DGetFlowVar (double*,double*,double*,double*,double*,double*,double*,void*);

double NavierStokes3DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims];
  double *u   = solver->u;

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double rho, vx, vy, vz, e, P, c, dxinv, dyinv, dzinv, local_cfl[3];
    IERR NavierStokes3DGetFlowVar(&u[nvars*p],&rho,&vx,&vy,&vz,&e,&P,param); CHECKERR(ierr);

    c         = sqrt(param->gamma*P/rho); /* speed of sound */
    dxinv     = solver->GetCoordinate(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv); /* 1/dx */
    dyinv     = solver->GetCoordinate(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv); /* 1/dy */
    dzinv     = solver->GetCoordinate(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->dxinv); /* 1/dz */

    local_cfl[_XDIR_] = (absolute(vx)+c)*dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vy)+c)*dt*dyinv; /* local cfl for this grid point (y) */
    local_cfl[_ZDIR_] = (absolute(vz)+c)*dt*dzinv; /* local cfl for this grid point (z) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];
    if (local_cfl[_ZDIR_] > max_cfl) max_cfl = local_cfl[_ZDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
