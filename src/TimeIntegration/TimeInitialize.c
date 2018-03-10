/*! @file TimeInitialize.c
    @brief Initialize time integration
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeRHSFunctionExplicit(double*,double*,void*,void*,double);

/*!
  Initialize time integration: This function initializes all that is required 
  for time integration. 
  + It sets the number of iterations, time step size, simulation time, etc.
  + It allocates solution, right-hand-side, and stage solution arrays needed
    by specific time integration methods.
  + It calls the method-specific initialization functions.
*/
int TimeInitialize(
                    void *s,  /*!< Solver object of type #HyPar */
                    void *m,  /*!< MPI object of type #MPIVariables */
                    void *ts  /*!< Time integration object of type #TimeIntegration */
                  )
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           s;
  MPIVariables    *mpi    = (MPIVariables*)    m;
  int             d, i;
  if (!solver) return(1);

  TS->solver        = solver;
  TS->mpi           = mpi;
  TS->n_iter        = solver->n_iter;
  TS->t_final       = solver->t_final;
  TS->restart_iter  = solver->restart_iter;
  TS->restart_time  = solver->restart_time;
  TS->dt            = solver->dt;
  TS->cfl           = solver->cfl;
  TS->max_cfl       = 0.0;
  TS->norm          = 0.0;
  TS->TimeIntegrate = solver->TimeIntegrate;

  TS->iter = TS->restart_iter;

  if ((TS->n_iter == -1) && (TS->t_final == -1)) {
    fprintf(stderr, "ERROR in TimeInitialize(): both n_iter and t_final cannot be unspecified.\n");
    return(1);
  }
  if ((TS->dt == -1) && (TS->cfl == -1)) {
    fprintf(stderr, "ERROR in TimeInitialize(): both dt and cfl cannot be unspecified.\n");
    return(1);
  }
  if ((TS->dt >= 0) && (TS->cfl >= 0)) {
    fprintf(stderr, "WARNING in TimeInitialize(): both dt and cfl have been specified. The simulation will continue with fixed CFL now.\n");
  }
  if ((TS->cfl >= 0) && (!solver->ComputeCFL)) {
    fprintf(stderr,"ERROR in TimeInitialize(): cfl specified but solver->ComputeCFL is NULL.\n");
    return(1);
  }

  if (TS->restart_time == 0) {
    TS->waqt = max(TS->restart_time, (double) TS->restart_iter * TS->dt);
  } else {
    TS->waqt = TS->restart_time;
  }

  int size = solver->nvars;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d] + 2*solver->ghosts);
  TS->u   = (double*) calloc (size,sizeof(double));
  TS->rhs = (double*) calloc (size,sizeof(double));
  _ArraySetValue_(TS->u  ,size,0.0);
  _ArraySetValue_(TS->rhs,size,0.0);

  /* initialize arrays to NULL, then allocate as necessary */
  TS->U             = NULL;
  TS->Udot          = NULL;
  TS->BoundaryFlux  = NULL;
  
  if (!strcmp(solver->time_scheme,_RK_)) {

    /* explicit Runge-Kutta methods */
    ExplicitRKParameters  *params = (ExplicitRKParameters*)  solver->msti;
    int nstages = params->nstages;
    TS->U     = (double**) calloc (nstages,sizeof(double*));
    TS->Udot  = (double**) calloc (nstages,sizeof(double*));
    for (i = 0; i < nstages; i++) {
      TS->U[i]    = (double*) calloc (size,sizeof(double));
      TS->Udot[i] = (double*) calloc (size,sizeof(double));
    }
    TS->BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
    for (i=0; i<nstages; i++) 
      TS->BoundaryFlux[i] = (double*) calloc (2*solver->ndims*solver->nvars,sizeof(double));
  
  } else if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
    
    /* General Linear Methods with Global Error Estimate */
    GLMGEEParameters *params = (GLMGEEParameters*) solver->msti;
    int nstages = params->nstages;
    int r       = params->r;
    TS->U     = (double**) calloc (2*r-1  ,sizeof(double*));
    TS->Udot  = (double**) calloc (nstages,sizeof(double*));
    for (i=0; i<2*r-1; i++)   TS->U[i]    = (double*) calloc (size,sizeof(double));
    for (i=0; i<nstages; i++) TS->Udot[i] = (double*) calloc (size,sizeof(double));
    TS->BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
    for (i=0; i<nstages; i++) 
      TS->BoundaryFlux[i] = (double*) calloc (2*solver->ndims*solver->nvars,sizeof(double));
    
    if (!strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
      for (i=0; i<r-1; i++) _ArrayCopy1D_(solver->u,TS->U[r+i],size);
    } else {
      for (i=0; i<r-1; i++) _ArraySetValue_(TS->U[r+i],size,0.0);
    }
  }

  /* set right-hand side function pointer */
  TS->RHSFunction = TimeRHSFunctionExplicit;

  /* open files for writing */
  if (!mpi->rank) {
    if (solver->write_residual) TS->ResidualFile = (void*) fopen("residual.out","w");
    else                        TS->ResidualFile = NULL;
  } else                        TS->ResidualFile = NULL;

  solver->time_integrator = TS;
  return(0);
}

