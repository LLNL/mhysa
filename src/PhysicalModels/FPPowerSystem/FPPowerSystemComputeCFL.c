#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <mpivars.h>
#include <hypar.h>

inline double FPPowerSystemDriftFunction(int,void*,double,double,double);

double FPPowerSystemComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar         *solver = (HyPar*)        s;
  FPPowerSystem *params = (FPPowerSystem*)solver->physics;
  int           ierr = 0;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_cfl = 0;
  int *index  = (int*) calloc (ndims,sizeof(int));
  int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    double x      = solver->GetCoordinate(0,index[0],dim,ghosts,solver->x);
    double y      = solver->GetCoordinate(1,index[1],dim,ghosts,solver->x);
    double dxinv  = solver->GetCoordinate(0,index[0],dim,ghosts,solver->dxinv);
    double dyinv  = solver->GetCoordinate(1,index[1],dim,ghosts,solver->dxinv);
    double drift_x= FPPowerSystemDriftFunction(0,params,x,y,t);
    double drift_y= FPPowerSystemDriftFunction(1,params,x,y,t);

    double local_cfl_x = absolute(drift_x) * dt * dxinv;
    double local_cfl_y = absolute(drift_y) * dt * dyinv;

    if (local_cfl_x > max_cfl) max_cfl = local_cfl_x;
    if (local_cfl_y > max_cfl) max_cfl = local_cfl_y;

    done = ArrayIncrementIndex(ndims,dim,index);
  }

  free(index);
  return(max_cfl);
}
