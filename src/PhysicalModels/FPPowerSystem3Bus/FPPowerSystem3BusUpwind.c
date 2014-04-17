#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

int FPPowerSystem3BusDriftFunction(int,void*,double*,double,double*);

int FPPowerSystem3BusUpwind(double *fI,double *fL,double *fR,double *uL,double *uR,
                        double *u,int dir,void *s,double t)
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->physics;
  int               done,v;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ghosts  = solver->ghosts;
  int *dim    = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0, p);
      double x[ndims]; /* coordinates of the interface */
      double x1, x2, drift[ndims];
      if (dir == 0) {
        _GetCoordinate_(0,index_inter[0]-1,dim,ghosts,solver->x,x1);
        _GetCoordinate_(0,index_inter[0]  ,dim,ghosts,solver->x,x2);
        x[0] = 0.5 * ( x1 + x2 );
        _GetCoordinate_(1,index_inter[1],dim,ghosts,solver->x,x[1]);
        _GetCoordinate_(2,index_inter[2],dim,ghosts,solver->x,x[2]);
        _GetCoordinate_(3,index_inter[3],dim,ghosts,solver->x,x[3]);
      } else if (dir == 1) {
        _GetCoordinate_(0,index_inter[0],dim,ghosts,solver->x,x[0]);
        _GetCoordinate_(1,index_inter[1]-1,dim,ghosts,solver->x,x1);
        _GetCoordinate_(1,index_inter[1]  ,dim,ghosts,solver->x,x2);
        x[1] = 0.5 * ( x1 + x2 );
        _GetCoordinate_(2,index_inter[2],dim,ghosts,solver->x,x[2]);
        _GetCoordinate_(3,index_inter[3],dim,ghosts,solver->x,x[3]);
      } else if (dir == 2) {
        _GetCoordinate_(0,index_inter[0],dim,ghosts,solver->x,x[0]);
        _GetCoordinate_(1,index_inter[1],dim,ghosts,solver->x,x[1]);
        _GetCoordinate_(2,index_inter[2]-1,dim,ghosts,solver->x,x1);
        _GetCoordinate_(2,index_inter[2]  ,dim,ghosts,solver->x,x2);
        x[2] = 0.5 * ( x1 + x2 );
        _GetCoordinate_(3,index_inter[3],dim,ghosts,solver->x,x[3]);
      } else if (dir == 3) {
        _GetCoordinate_(0,index_inter[0],dim,ghosts,solver->x,x[0]);
        _GetCoordinate_(1,index_inter[1],dim,ghosts,solver->x,x[1]);
        _GetCoordinate_(2,index_inter[2],dim,ghosts,solver->x,x[2]);
        _GetCoordinate_(3,index_inter[3]-1,dim,ghosts,solver->x,x1);
        _GetCoordinate_(3,index_inter[3]  ,dim,ghosts,solver->x,x2);
        x[3] = 0.5 * ( x1 + x2 );
      }
      FPPowerSystem3BusDriftFunction(dir,params,x,t,drift);
      for (v = 0; v < nvars; v++)  
        fI[nvars*p+v] = drift[dir] * (drift[dir] > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
