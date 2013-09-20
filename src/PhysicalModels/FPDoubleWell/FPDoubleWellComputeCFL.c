#include <mathfunctions.h>
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
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;

  int     offset  = 0;
  double  max_cfl = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double advection_speed = 0;
        double x = solver->x[i+ghosts];
        double local_cfl =  absolute(drift(x)) * dt 
                          * dxinv[offset+ghosts+i];
        if (local_cfl > max_cfl) max_cfl = local_cfl;
      }
    }
    offset += dim[d] + 2*ghosts;
  }

  return(max_cfl);
}
