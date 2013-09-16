#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeForwardEuler(void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;
  int             ierr    = 0, d;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,solver->u);     CHECKERR(ierr);
  ierr = MPIExchangeBoundaries  (solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,solver->u);     CHECKERR(ierr);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  double *rhs = TS->rhs;
  ierr = ArraySetValue_double(rhs,size*solver->nvars,0.0);          CHECKERR(ierr);
  if (solver->HyperbolicFunction) {
    ierr = solver->HyperbolicFunction(solver,mpi);                  CHECKERR(ierr);
    ierr = ArrayAXPY(solver->hyp    ,-1.0,rhs,size*solver->nvars);  CHECKERR(ierr);
  }
  if (solver->ParabolicFunction) {
    ierr = solver->ParabolicFunction (solver,mpi);                  CHECKERR(ierr);
    ierr = ArrayAXPY(solver->par    , 1.0,rhs,size*solver->nvars);  CHECKERR(ierr);
  }
  if (solver->SourceFunction) {
    ierr = solver->SourceFunction    (solver,mpi);                  CHECKERR(ierr);
    ierr = ArrayAXPY(solver->source , 1.0,rhs,size*solver->nvars);  CHECKERR(ierr);
  }

  /* Update solution */
  ierr = ArrayAXPY(rhs,solver->dt,solver->u,size*solver->nvars);    CHECKERR(ierr);

  return(0);
}
