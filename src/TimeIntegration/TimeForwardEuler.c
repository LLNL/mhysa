/*! @file TimeForwardEuler.c
    @brief Forward Euler method
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Advance the ODE given by
  \f{equation}{
    \frac{d{\bf u}}{dt} = {\bf F} \left({\bf u}\right)
  \f}
  by one time step of size #TimeIntegration::dt using the forward Euler method
  given by
  \f{equation}{
    {\bf u}^{n+1} = {\bf u}^n + \Delta t {\bf F}\left( {\bf u}^n \right)
  \f}
  where the superscript represents the time level, \f$\Delta t\f$ is the
  time step size #TimeIntegration::dt, and \f${\bf F}\left({\bf u}\right)\f$ is 
  computed by #TimeIntegration::RHSFunction.
*/
int TimeForwardEuler(
                      void *ts /*!< Time integrator object of type #TimeIntegration */
                    )
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
  _ArrayAXPY_(TS->rhs,TS->dt,solver->u,size*solver->nvars);
  _ArrayScaleCopy1D_(solver->StageBoundaryIntegral,TS->dt,
                     solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars);

  if (solver->PostStage) 
    { IERR solver->PostStage(solver->u,solver,mpi,TS->waqt);     CHECKERR(ierr); }

  return(0);
}
