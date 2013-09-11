#include <stdio.h>
#include <stdlib.h>
#include <mpivars.h>
#include <hypar.h>

int Cleanup(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0;

  if (!mpi->rank) printf("Deallocating arrays.\n");

  /* These variables are allocated in Initialize.c */
  free(solver->dim_global);
  free(solver->dim_local);
  free(solver->index);
  free(solver->u);
  free(solver->x);
  free(mpi->iproc);
  free(mpi->ip);
  free(mpi->is);
  free(mpi->ie);

  return(0);
}
