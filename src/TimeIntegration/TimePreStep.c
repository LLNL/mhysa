/*! @file TimePreStep.c
    @brief Pre-time-step function
    @author Debojyoti Ghosh
*/

#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*! 
  Pre-time-step function: This function is called before each time
  step. Some notable things this does are:
  + Computes CFL and diffusion numbers.
  + Call the physics-specific pre-time-step function, if defined.
*/
int TimePreStep(void *ts /*!< Object of type #TimeIntegration */ )
{
  TimeIntegration *TS      = (TimeIntegration*) ts;
  HyPar           *solver  = (HyPar*)           TS->solver;
  MPIVariables    *mpi     = (MPIVariables*)    TS->mpi;
  _DECLARE_IERR_;

  gettimeofday(&TS->iter_start_time,NULL);

  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,solver->u,NULL,TS->waqt);   CHECKERR(ierr);
  IERR solver->ApplyIBConditions(solver,mpi,solver->u,TS->waqt);              CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,solver->u);               CHECKERR(ierr);

  int size = 1,d;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d] + 2*solver->ghosts);
  _ArrayCopy1D_(solver->u,TS->u,size*solver->nvars);

  /* compute max stable dt over the domain */
  double local_max_cfl = -1.0, global_max_cfl = -1.0;
  if (solver->ComputeCFL) local_max_cfl = solver->ComputeCFL(solver,mpi,1.0,TS->waqt);
  IERR MPIMax_double(&global_max_cfl,&local_max_cfl,1,&mpi->world); CHECKERR(ierr);
  double stable_dt = 1.0/global_max_cfl;

  /* if CFL-based time stepping, then compute dt based on CFL */
  if (TS->cfl >= 0) {
    TS->dt = TS->cfl * stable_dt;
  }
  /* if final step, adjust dt to match final simulation time */
  if ((TS->t_final >= 0) && ((TS->waqt+TS->dt) > TS->t_final)) {
    TS->dt = TS->t_final - TS->waqt;
  }
  /* compute CFL of current step */
  TS->max_cfl = TS->dt / stable_dt;

  /* set the step boundary flux integral value to zero */
  _ArraySetValue_(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);

  if (solver->PreStep)  { IERR solver->PreStep(solver->u,solver,mpi,TS->waqt); CHECKERR(ierr); }

  return(0);
}

