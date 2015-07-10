#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

int FPPowerSystem1BusPostStep(double *u,void* s,void *m,double t,int iter)
{
  HyPar             *solver = (HyPar*)              s;
  MPIVariables      *mpi    = (MPIVariables*)       m;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*)  solver->physics;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int index[ndims];

  double local_sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dxinv; _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->dxinv,dxinv);
    double dyinv; _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->dxinv,dyinv);
    local_sum     += (u[p] / (dxinv * dyinv));
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  double local_integral = local_sum;
  double global_integral = 0; 
  IERR MPISum_double(&global_integral,&local_integral,1,&mpi->world); CHECKERR(ierr);
  params->pdf_integral = global_integral;

  return(0);
}
