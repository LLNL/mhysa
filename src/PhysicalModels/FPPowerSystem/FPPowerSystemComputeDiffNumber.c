#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <mpivars.h>
#include <hypar.h>

inline double FPPowerSystemDissipationFunction(int,void*,double);

double FPPowerSystemComputeDiffNumber(void *s,void *m,double dt,double t)
{
  HyPar         *solver = (HyPar*)        s;
  FPPowerSystem *params = (FPPowerSystem*)solver->physics;
  int           ierr = 0;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_diff = 0;
  int     index[ndims];
  int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    double dxinv  = solver->GetCoordinate(0,index[0],dim,ghosts,solver->dxinv);
    double dyinv  = solver->GetCoordinate(1,index[1],dim,ghosts,solver->dxinv);
    double dissp_x= FPPowerSystemDissipationFunction(0,params,t);
    double dissp_y= FPPowerSystemDissipationFunction(1,params,t);

    double local_diff_x = absolute(dissp_x) * dt * dxinv * dxinv;
    double local_diff_y = absolute(dissp_y) * dt * dyinv * dyinv;

    if (local_diff_x > max_diff) max_diff = local_diff_x;
    if (local_diff_y > max_diff) max_diff = local_diff_y;

    done = ArrayIncrementIndex(ndims,dim,index);
  }

  return(max_diff);
}
