#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystem1BusDissipationFunction(int,int,void*,double);

double FPPowerSystem1BusComputeDiffNumber(void *s,void *m,double dt,double t)
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*)  solver->physics;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_diff = 0;
  int     index[ndims];
  int     done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv);
    double dyinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dyinv);
    double dissp_yx= FPPowerSystem1BusDissipationFunction(_YDIR_,_XDIR_,params,t);
    double dissp_yy= FPPowerSystem1BusDissipationFunction(_YDIR_,_YDIR_,params,t);

    double local_diff_yx = absolute(dissp_yx) * dt * dyinv * dxinv;
    double local_diff_yy = absolute(dissp_yy) * dt * dyinv * dyinv;

    if (local_diff_yx > max_diff) max_diff = local_diff_yx;
    if (local_diff_yy > max_diff) max_diff = local_diff_yy;

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_diff);
}
