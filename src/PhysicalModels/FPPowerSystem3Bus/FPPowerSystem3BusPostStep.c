#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

int FPPowerSystem3BusPostStep(double *u,void* s,void *m,double t)
{
  HyPar             *solver = (HyPar*)              s;
  MPIVariables      *mpi    = (MPIVariables*)       m;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->physics;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int index[ndims];

  double local_sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dxinv[ndims];
    _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv[0]);
    _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dxinv[1]);
    _GetCoordinate_(2,index[2],dim,ghosts,solver->dxinv,dxinv[2]);
    _GetCoordinate_(3,index[3],dim,ghosts,solver->dxinv,dxinv[3]);
    local_sum     += (u[p] / (dxinv[0]*dxinv[1]*dxinv[2]*dxinv[3]));
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  double local_integral = local_sum;
  double global_integral = 0; 
  IERR MPISum_double(&global_integral,&local_integral,1,&mpi->world); CHECKERR(ierr);
  params->pdf_integral = global_integral;

  return(0);
}
