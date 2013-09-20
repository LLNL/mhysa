#include <fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

double FPDoubleWellComputeCFL(void *s,void *m,double dt)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  FPDoubleWell  *params = (FPDoubleWell*) solver->physics;
  int           ierr    = 0, d, i, v;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;

  int     offset  = 0;
  double  max_cfl = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double advection_speed = 0;
        double x = solver->x[i];
        advection_speed = 4*x*(x*x-1.0);
        double local_cfl = advection_speed*dt*dxinv[offset+i];
        if (local_cfl > max_cfl) max_cfl = local_cfl;
      }
    }
    offset += dim[d];
  }

  return(max_cfl);
}
