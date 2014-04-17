#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <mpivars.h>
#include <hypar.h>

int FPPowerSystem3BusDriftFunction(int,void*,double*,double,double*);

double FPPowerSystem3BusComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->physics;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_cfl = 0;
  int     index[ndims];
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    double x[ndims], dxinv[ndims], drift[ndims];
    _GetCoordinate_(0,index[0],dim,ghosts,solver->x,x[0]);
    _GetCoordinate_(1,index[1],dim,ghosts,solver->x,x[1]);
    _GetCoordinate_(2,index[2],dim,ghosts,solver->x,x[2]);
    _GetCoordinate_(3,index[3],dim,ghosts,solver->x,x[3]);
    _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv[0]);
    _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dxinv[1]);
    _GetCoordinate_(2,index[2],dim,ghosts,solver->dxinv,dxinv[2]);
    _GetCoordinate_(3,index[3],dim,ghosts,solver->dxinv,dxinv[3]);

    FPPowerSystem3BusDriftFunction(0,params,x,t,drift);

    double local_cfl[ndims];
    local_cfl[0] = absolute(drift[0]) * dt * dxinv[0];
    local_cfl[1] = absolute(drift[1]) * dt * dxinv[1];
    local_cfl[2] = absolute(drift[2]) * dt * dxinv[2];
    local_cfl[3] = absolute(drift[3]) * dt * dxinv[3];

    if (local_cfl[0] > max_cfl) max_cfl = local_cfl[0];
    if (local_cfl[1] > max_cfl) max_cfl = local_cfl[1];
    if (local_cfl[2] > max_cfl) max_cfl = local_cfl[2];
    if (local_cfl[3] > max_cfl) max_cfl = local_cfl[3];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
