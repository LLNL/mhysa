#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimePreStep(void *ts)
{
  TimeIntegration *TS      = (TimeIntegration*) ts;
  HyPar           *solver  = (HyPar*)           TS->solver;
  MPIVariables    *mpi     = (MPIVariables*)    TS->mpi;
  int             ierr     = 0;

  /* copy current solution for norm computation later */
  int size = 1,d;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d] + 2*solver->ghosts);
  ierr = ArrayCopy1D_double(solver->u,TS->u,size*solver->nvars); CHECKERR(ierr);

  /* compute max CFL over the domain */
  double local_max_cfl = -1.0;
  if (solver->ComputeCFL) local_max_cfl = solver->ComputeCFL(solver,mpi,TS->dt);
  ierr = MPIMax_double(&TS->max_cfl,&local_max_cfl,1); CHECKERR(ierr);

  return(0);
}

