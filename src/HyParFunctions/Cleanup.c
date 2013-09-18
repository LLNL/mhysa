#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <timeintegration.h>
#include <hypar.h>
#include <physics.h>

/* include header files for each physical model */
#include <advectiondiffusionreaction.h>

int Cleanup(void *s,void *m)
{
  HyPar           *solver   = (HyPar*)          s;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  int             ierr    = 0,i;

  if (!mpi->rank) printf("Deallocating arrays.\n");

  /* Clean up boundary zones */
  for (i = 0; i < solver->nBoundaryZones; i++) {
    ierr = BCCleanup(&boundary[i]); CHECKERR(ierr);
  }
  free(solver->boundary);

  /* Clean up any allocations in physical model */
  if (!strcmp(solver->model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {
    ierr = LinearADRCleanup(solver->physics); CHECKERR(ierr);
  }
  free(solver->physics);

  /* Clean up any allocations from time-integration */
  if (solver->msti) {
    ierr = TimeMSTICleanup(solver->msti); CHECKERR(ierr);
    free(solver->msti);
  }

  /* Clean up any spatial reconstruction related allocations */
  if (solver->interp) free(solver->interp);

  /* These variables are allocated in Initialize.c */
  free(solver->dim_global);
  free(solver->dim_local);
  free(solver->index);
  free(solver->u);
  free(solver->hyp);
  free(solver->par);
  free(solver->source);
  free(solver->x);
  free(solver->dxinv);
  free(mpi->iproc);
  free(mpi->ip);
  free(mpi->is);
  free(mpi->ie);

  return(0);
}
