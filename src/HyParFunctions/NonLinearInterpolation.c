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

  int flag = (FluxFunction && solver->flag_nonlinearinterp && solver->SetInterpLimiterVar);
  if (flag) {;
    int     ndims  = solver->ndims;
    int     nvars  = solver->nvars;
    int     ghosts = solver->ghosts;
    int     *dim   = solver->dim_local;
    double  *x     = solver->x;

    int size = 1, d;
    for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

    /* apply boundary conditions and exchange data over MPI interfaces */
    IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);      CHECKERR(ierr);
    IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u);       CHECKERR(ierr);

    int offset = 0;
    for (d = 0; d < ndims; d++) {
      /* evaluate cell-centered flux */
      IERR FluxFunction(FluxC,u,d,solver,t);                          CHECKERR(ierr);
      /* calculate non-linear interpolation coefficients */
      IERR solver->SetInterpLimiterVar(FluxC,u,x+offset,d,solver,mpi);CHECKERR(ierr);
      offset += dim[d] + 2*ghosts;
    }
  }

  return(0);
}

