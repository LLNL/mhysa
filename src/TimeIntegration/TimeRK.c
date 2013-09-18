#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

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
    ierr = TS->RHSFunction(TS->Udot[stage],TS->U[stage],solver,mpi);
  }

  /* Stage completion */
  for (stage = 0; stage < params->nstages; stage++) {
    ierr = ArrayAXPY(TS->Udot[stage],solver->dt*params->b[stage],solver->u,size*solver->nvars);
    CHECKERR(ierr);
  }

  return(0);
}

