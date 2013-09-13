#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>

int Cleanup(void *s,void *m)
{
  HyPar           *solver   = (HyPar*)          s;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  int           ierr    = 0,i;

  if (!mpi->rank) printf("Deallocating arrays.\n");

  /* Clean up boundary zones */
  for (i = 0; i < solver->nBoundaryZones; i++) {
    ierr = BCCleanup(&boundary[i]); CHECKERR(ierr);
  }

  /* These variables are allocated in Initialize.c */
  free(solver->dim_global);
  free(solver->dim_local);
  free(solver->index);
  free(solver->u);
  free(solver->hyp);
  free(solver->par);
  free(solver->x);
  free(mpi->iproc);
  free(mpi->ip);
  free(mpi->is);
  free(mpi->ie);

  return(0);
}
