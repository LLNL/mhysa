#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

int FPPowerSystemPostStep(double *u,void* s,void *m,double t)
{
  HyPar         *solver = (HyPar*)         s;
  MPIVariables  *mpi    = (MPIVariables*)  m;
  FPPowerSystem *params = (FPPowerSystem*) solver->physics;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int index[ndims];

  double local_sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv);
    double dyinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dyinv);
    local_sum     += (u[p] / (dxinv * dyinv));
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  double local_integral = local_sum;
  double global_integral = 0; 
  IERR MPISum_double(&global_integral,&local_integral,1,&mpi->world); CHECKERR(ierr);
  params->pdf_integral = global_integral;

  return(0);
}
