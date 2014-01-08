#include <stdlib.h>
#include <basic.h>
#include <math.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

double Euler1DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int index[ndims];
  double *u   = solver->u;

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double rho, v, e, P, c, dxinv, local_cfl;
    _Euler1DGetFlowVar_((u+_MODEL_NVARS_*p),rho,v,e,P,param);
    _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */
    c         = sqrt(param->gamma*P/rho); /* speed of sound */
    local_cfl = (absolute(v)+c)*dt*dxinv; /* local cfl for this grid point */
    if (local_cfl > max_cfl) max_cfl = local_cfl;
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
