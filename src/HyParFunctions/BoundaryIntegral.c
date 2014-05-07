#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int BoundaryIntegral(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int d,v;

  double *local_integral  = (double*) calloc (nvars,sizeof(double));
  double *global_integral = (double*) calloc (nvars,sizeof(double));

  /* calculate the local boundary integral on each process */
  _ArraySetValue_(local_integral,nvars,0.0);
  for (d=0; d<ndims; d++) {
    for (v=0; v<nvars; v++) {
      local_integral[v] += (solver->StepBoundaryIntegral[(2*d+0)*nvars+v]);
      local_integral[v] += (solver->StepBoundaryIntegral[(2*d+1)*nvars+v]);
    }
  }

  /* add across process to calculate global boundary integral 
   * (internal (MPI) boundaries must cancel out                 */
  IERR MPISum_double(global_integral,local_integral,nvars,&mpi->world); CHECKERR(ierr);

  /* add to the total boundary integral */
  _ArrayAXPY_(global_integral,1.0,solver->TotalBoundaryIntegral,nvars);

  free(local_integral);
  free(global_integral);
  return(0);
}
