/*! @file TimeRK.c
    @brief Explicit Runge-Kutta method
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
  \f{align}{
    {\bf U}^{\left(i\right)} &= {\bf u}_n + \Delta t \sum_{j=1}^{i-1} a_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), \\
    {\bf u}_{n+1} &= {\bf u}_n + \Delta t \sum_{i=1}^s b_{i} {\bf F}\left({\bf U}^{\left(i\right)}\right),
  \f}
  where the subscript represents the time level, the superscripts represent the stages, \f$\Delta t\f$ is the
  time step size #TimeIntegration::dt, and \f${\bf F}\left({\bf u}\right)\f$ is computed by #TimeIntegration::RHSFunction.
  The Butcher tableaux coefficients are \f$a_{ij}\f$ (#ExplicitRKParameters::A) and \f$b_i\f$ 
  (#ExplicitRKParameters::b).

  Note: In the code #TimeIntegration::Udot is equivalent to \f${\bf F}\left({\bf u}\right)\f$.
*/
int TimeRK(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration       *TS     = (TimeIntegration*) ts;
  HyPar                 *solver = (HyPar*)           TS->solver;
  MPIVariables          *mpi    = (MPIVariables*)    TS->mpi;
  ExplicitRKParameters  *params = (ExplicitRKParameters*)  solver->msti;
  int                   d, stage, i;
  _DECLARE_IERR_;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* Calculate stage values */
  for (stage = 0; stage < params->nstages; stage++) {
    double stagetime = TS->waqt + params->c[stage]*TS->dt;
    _ArrayCopy1D_(solver->u,TS->U[stage],size*solver->nvars);
    for (i = 0; i < stage; i++) {
      _ArrayAXPY_(TS->Udot[i],TS->dt*params->A[stage*params->nstages+i],
                  TS->U[stage],size*solver->nvars); 
    }
    if (solver->PreStage)
      { IERR solver->PreStage(stage,TS->U ,solver,mpi,stagetime); CHECKERR(ierr); }
    IERR TS->RHSFunction(TS->Udot[stage],TS->U[stage],solver,mpi,stagetime);
    if (solver->PostStage) 
      { IERR solver->PostStage(TS->U[stage],solver,mpi,stagetime); CHECKERR(ierr); }

    _ArraySetValue_(TS->BoundaryFlux[stage],2*solver->ndims*solver->nvars,0.0);
    _ArrayCopy1D_(solver->StageBoundaryIntegral,TS->BoundaryFlux[stage],2*solver->ndims*solver->nvars);
  }

  /* Step completion */
  for (stage = 0; stage < params->nstages; stage++) {
    _ArrayAXPY_(TS->Udot[stage],TS->dt*params->b[stage],solver->u,size*solver->nvars);
    _ArrayAXPY_(TS->BoundaryFlux[stage],TS->dt*params->b[stage],solver->StepBoundaryIntegral,
                2*solver->ndims*solver->nvars);
  }

  return(0);
}

