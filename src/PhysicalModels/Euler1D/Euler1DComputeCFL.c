#include <stdlib.h>
#include <math.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;

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
    double rho, v, e, P, c, dxinv, local_cfl;
    _Euler1DGetFlowVar_((u+nvars*p),rho,v,e,P,param);

    c         = sqrt(param->gamma*P/rho); /* speed of sound */
    dxinv     = solver->GetCoordinate(0,index[0],dim,ghosts,solver->dxinv); /* 1/dx */
    local_cfl = (absolute(v)+c)*dt*dxinv; /* local cfl for this grid point */
    if (local_cfl > max_cfl) max_cfl = local_cfl;

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(0);
}
