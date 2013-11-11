#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

inline int Euler1DGetFlowVar (double*,double*,double*,double*,double*,void*);

int Euler1DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;
  int       ierr    = 0;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims];
  double *u   = solver->u;

  double max_cfl = 0;
  int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    double rho, v, e, P, c, dxinv, local_cfl;
    ierr = Euler1DGetFlowVar(&u[nvars*p],&rho,&v,&e,&P,param); CHECKERR(ierr);

    c         = sqrt(param->gamma*P/rho); /* speed of sound */
    dxinv     = solver->GetCoordinate(0,index[0],dim,ghosts,solver->dxinv); /* 1/dx */
    local_cfl = (absolute(v)+c)*dt*dxinv; /* local cfl for this grid point */
    if (local_cfl > max_cfl) max_cfl = local_cfl;

    done = ArrayIncrementIndex(ndims,dim,index);
  }

  return(0);
}
