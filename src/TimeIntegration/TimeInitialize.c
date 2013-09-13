#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeInitialize(void *s,void *m,void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           s;
  MPIVariables    *mpi    = (MPIVariables*)    m;
  int             ierr    = 0, d;
  if (!solver) return(1);

  TS->solver = solver;
  TS->mpi    = mpi;
  TS->n_iter = solver->n_iter;
  TS->waqt   = 0.0;
  TS->dt     = solver->dt;
  TS->max_cfl= 0.0;
  TS->norm   = 0.0;
  TS->TimeIntegrate = solver->TimeIntegrate;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d] + 2*solver->ghosts);
  TS->u = (double*) calloc (size*solver->nvars,sizeof(double));
  ierr = ArraySetValue_double(TS->u,size*solver->nvars,0.0); CHECKERR(ierr);

  /* open files for writing */
  if (!mpi->rank) {
    if (solver->write_residual) TS->ResidualFile = (void*) fopen("residual.out","w");
    else                        TS->ResidualFile = NULL;
  } else                        TS->ResidualFile = NULL;

  return(0);
}

