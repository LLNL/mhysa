#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeRHSFunctionExplicit(double*,double*,void*,void*,double);

int TimeInitialize(void *s,void *m,void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           s;
  MPIVariables    *mpi    = (MPIVariables*)    m;
  MSTIParameters  *params = (MSTIParameters*)  solver->msti;
  int             d;
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
  TS->u   = (double*) calloc (size*solver->nvars,sizeof(double));
  TS->rhs = (double*) calloc (size*solver->nvars,sizeof(double));
  _ArraySetValue_(TS->u  ,size*solver->nvars,0.0);
  _ArraySetValue_(TS->rhs,size*solver->nvars,0.0);

  /* allocate arrays for multi-stage schemes, if required */
  if (params) {
    int nstages = params->nstages;
    TS->U     = (double**) calloc (nstages,sizeof(double*));
    TS->Udot  = (double**) calloc (nstages,sizeof(double*));
    int i;
    for (i = 0; i < nstages; i++) {
      TS->U[i]    = (double*) calloc (size*solver->nvars,sizeof(double));
      TS->Udot[i] = (double*) calloc (size*solver->nvars,sizeof(double));
    }
  } else TS->U = TS->Udot = NULL;

  /* set right-hand side function pointer */
  TS->RHSFunction = TimeRHSFunctionExplicit;

  /* open files for writing */
  if (!mpi->rank) {
    if (solver->write_residual) TS->ResidualFile = (void*) fopen("residual.out","w");
    else                        TS->ResidualFile = NULL;
  } else                        TS->ResidualFile = NULL;

  return(0);
}

