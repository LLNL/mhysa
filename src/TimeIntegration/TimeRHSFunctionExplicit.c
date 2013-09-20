#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int TimeRHSFunctionExplicit(double *rhs,double *u,void *s,void *m, double t) 
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             ierr    = 0, d;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,u);               CHECKERR(ierr);
  ierr = MPIExchangeBoundaries  (solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,u);               CHECKERR(ierr);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  ierr = ArraySetValue_double(rhs,size*solver->nvars,0.0);            CHECKERR(ierr);
  if (solver->HyperbolicFunction) {
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t);    CHECKERR(ierr);
    ierr = ArrayAXPY(solver->hyp    ,-1.0,rhs,size*solver->nvars);    CHECKERR(ierr);
  }
  if (solver->ParabolicFunction) {
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);    CHECKERR(ierr);
    ierr = ArrayAXPY(solver->par    , 1.0,rhs,size*solver->nvars);    CHECKERR(ierr);
  }
  if (solver->SourceFunction) {
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t); CHECKERR(ierr);
    ierr = ArrayAXPY(solver->source , 1.0,rhs,size*solver->nvars);    CHECKERR(ierr);
  }

  return(0);
}
