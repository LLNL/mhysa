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
  int             d;
  _DECLARE_IERR_;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  if (solver->PreStage) 
    { IERR solver->PreStage(0,&solver->u,solver,mpi,TS->waqt);   CHECKERR(ierr); }

  /* Evaluate right-hand side and update solution */
  IERR TS->RHSFunction(TS->rhs,solver->u,solver,mpi,TS->waqt);   CHECKERR(ierr);
  _ArrayAXPY_(TS->rhs,solver->dt,solver->u,size*solver->nvars);

  if (solver->PostStage) 
    { IERR solver->PostStage(solver->u,solver,mpi,TS->waqt);     CHECKERR(ierr); }

  return(0);
}
