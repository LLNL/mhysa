#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int NonLinearInterpolation(double *u,void *s,void *m,double t,int(*FluxFunction)(double*,double*,int,void*,double))
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  double        *FluxC  = solver->fluxC; /* cell centered flux */
  _DECLARE_IERR_;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *x     = solver->x;

  int size = 1, d;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  if (!FluxFunction) return(0);

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    /* evaluate cell-centered flux */
    IERR FluxFunction(FluxC,u,d,solver,t); CHECKERR(ierr);
    /* calculate non-linear interpolation coefficients */
    IERR solver->SetInterpLimiterVar(FluxC,u,x+offset,d,solver,mpi);
    offset += dim[d] + 2*ghosts;
  }

  return(0);
}

