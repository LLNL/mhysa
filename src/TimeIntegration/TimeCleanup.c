#include <stdio.h>
#include <stdlib.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeCleanup(void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;
  MSTIParameters  *params = (MSTIParameters*)  solver->msti;

  /* close files opened for writing */
  if (!mpi->rank) if (solver->write_residual) fclose((FILE*)TS->ResidualFile);

  /* deallocate arrays for multi-stage schemes */
  if (TS->U) {
    int i; for (i=0; i<params->nstages; i++) free(TS->U[i]);
    free(TS->U);
  }
  if (TS->Udot) {
    int i; for (i=0; i<params->nstages; i++) free(TS->Udot[i]);
    free(TS->Udot);
  }

  if (TS->BoundaryFlux) {
    int i; for (i=0; i<params->nstages; i++) free(TS->BoundaryFlux[i]);
    free(TS->BoundaryFlux);
  }

  /* deallocate arrays */
  free(TS->u  );
  free(TS->rhs);
  return(0);
}
