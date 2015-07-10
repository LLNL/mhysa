#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <physicalmodels/fpdoublewell.h>
#include <hypar.h>

int FPDoubleWellPostStep(double *u,void* s,void *m,double t,int iter)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  FPDoubleWell  *params = (FPDoubleWell*) solver->physics;
  _DECLARE_IERR_;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int index[ndims];

  double local_sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dx = 1.0 / solver->dxinv[index[0]+ghosts];
    local_sum     += (u[p] * dx);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  double local_integral = local_sum;
  double global_integral = 0; 
  IERR MPISum_double(&global_integral,&local_integral,1,&mpi->world); CHECKERR(ierr);
  params->pdf_integral = global_integral;

  return(0);
}
