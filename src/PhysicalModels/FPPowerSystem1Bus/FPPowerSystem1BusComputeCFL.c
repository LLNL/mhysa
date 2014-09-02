#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystem1BusDriftFunction(int,void*,double,double,double);

double FPPowerSystem1BusComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar         *solver = (HyPar*)        s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*)solver->physics;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_cfl = 0;
  int     index[ndims];
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    double x;     _GetCoordinate_(0,index[0],dim,ghosts,solver->x,x);
    double y;     _GetCoordinate_(1,index[1],dim,ghosts,solver->x,y);
    double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv);
    double dyinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dyinv);
    double drift_x= FPPowerSystem1BusDriftFunction(0,params,x,y,t);
    double drift_y= FPPowerSystem1BusDriftFunction(1,params,x,y,t);

    double local_cfl_x = absolute(drift_x) * dt * dxinv;
    double local_cfl_y = absolute(drift_y) * dt * dyinv;

    if (local_cfl_x > max_cfl) max_cfl = local_cfl_x;
    if (local_cfl_y > max_cfl) max_cfl = local_cfl_y;

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
