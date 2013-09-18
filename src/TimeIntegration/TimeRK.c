#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

static int ComputeRHS(double*,double*,void*,void*);

int TimeRK(void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;
  MSTIParameters  *params = (MSTIParameters*)  solver->msti;
  int             ierr    = 0, d, stage, i;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* Calculate stage values */
  for (stage = 0; stage < params->nstages; stage++) {
    ierr = ArrayCopy1D_double(solver->u,TS->U[stage],size*solver->nvars); CHECKERR(ierr);
    for (i = 0; i < stage; i++) {
      ierr = ArrayAXPY(TS->Udot[i],solver->dt*params->A[stage*params->nstages+i],
                       TS->U[stage],size*solver->nvars); 
      CHECKERR(ierr);
    }
    ierr = ComputeRHS(TS->Udot[stage],TS->U[stage],solver,mpi);
  }

  /* Stage completion */
  for (stage = 0; stage < params->nstages; stage++) {
    ierr = ArrayAXPY(TS->Udot[stage],solver->dt*params->b[stage],solver->u,size*solver->nvars);
    CHECKERR(ierr);
  }

  return(0);
}

int ComputeRHS(double *rhs,double *u,void *s,void *m) 
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             ierr    = 0, d;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,u);             CHECKERR(ierr);
  ierr = MPIExchangeBoundaries  (solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,u);             CHECKERR(ierr);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  ierr = ArraySetValue_double(rhs,size*solver->nvars,0.0);          CHECKERR(ierr);
  if (solver->HyperbolicFunction) {
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi);    CHECKERR(ierr);
    ierr = ArrayAXPY(solver->hyp    ,-1.0,rhs,size*solver->nvars);  CHECKERR(ierr);
  }
  if (solver->ParabolicFunction) {
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi);    CHECKERR(ierr);
    ierr = ArrayAXPY(solver->par    , 1.0,rhs,size*solver->nvars);  CHECKERR(ierr);
  }
  if (solver->SourceFunction) {
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi); CHECKERR(ierr);
    ierr = ArrayAXPY(solver->source , 1.0,rhs,size*solver->nvars);  CHECKERR(ierr);
  }

  return(0);
}
