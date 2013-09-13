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

  /* close files opened for writing */
  if (!mpi->rank) if (solver->write_residual) fclose((FILE*)TS->ResidualFile);

  /* deallocate arrays */
  free(TS->u);
  return(0);
}
