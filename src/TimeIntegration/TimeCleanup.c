/*! @file TimeCleanup.c
    @brief Clean up time integration
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Clean up all allocations related to time integration 
*/
int TimeCleanup(void *ts /*!< Object of type #TimeIntegration*/)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;

  /* close files opened for writing */
  if (!mpi->rank) if (solver->write_residual) fclose((FILE*)TS->ResidualFile);

  if (!strcmp(solver->time_scheme,_RK_)) {
    int i;
    ExplicitRKParameters  *params = (ExplicitRKParameters*)  solver->msti;
    for (i=0; i<params->nstages; i++) free(TS->U[i]);            free(TS->U);
    for (i=0; i<params->nstages; i++) free(TS->Udot[i]);         free(TS->Udot);
    for (i=0; i<params->nstages; i++) free(TS->BoundaryFlux[i]); free(TS->BoundaryFlux);
  } else if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
    int i;
    GLMGEEParameters  *params = (GLMGEEParameters*)  solver->msti;
    for (i=0; i<2*params->r-1  ; i++) free(TS->U[i]);            free(TS->U);
    for (i=0; i<params->nstages; i++) free(TS->Udot[i]);         free(TS->Udot);
    for (i=0; i<params->nstages; i++) free(TS->BoundaryFlux[i]); free(TS->BoundaryFlux);
  }

  /* deallocate arrays */
  free(TS->u  );
  free(TS->rhs);
  solver->time_integrator = NULL;
  return(0);
}
