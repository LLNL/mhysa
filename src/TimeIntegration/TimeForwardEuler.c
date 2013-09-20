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

  /* Evaluate right-hand side and update solution */
  ierr = TS->RHSFunction(TS->rhs,solver->u,solver,mpi,TS->waqt);            CHECKERR(ierr);
  ierr = ArrayAXPY(TS->rhs,solver->dt,solver->u,size*solver->nvars);        CHECKERR(ierr);

  return(0);
}
