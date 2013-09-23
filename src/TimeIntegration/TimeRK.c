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

  if (solver->PreStep)  { ierr = solver->PreStep(solver->u,solver,mpi,TS->waqt); CHECKERR(ierr); }

  /* Calculate stage values */
  for (stage = 0; stage < params->nstages; stage++) {
    double stagetime = TS->waqt + params->c[stage]*TS->dt;
    if (solver->PreStage) { 
      if (stage) { ierr = solver->PreStage(stage,TS->U ,solver,mpi,stagetime); CHECKERR(ierr); }
      else       { ierr = solver->PreStage(1,&solver->u,solver,mpi,stagetime); CHECKERR(ierr); }
    }
    ierr = ArrayCopy1D_double(solver->u,TS->U[stage],size*solver->nvars); CHECKERR(ierr);
    for (i = 0; i < stage; i++) {
      ierr = ArrayAXPY(TS->Udot[i],solver->dt*params->A[stage*params->nstages+i],
                       TS->U[stage],size*solver->nvars); 
      CHECKERR(ierr);
    }
    ierr = TS->RHSFunction(TS->Udot[stage],TS->U[stage],solver,mpi,stagetime);
    if (solver->PostStage) 
      { ierr = solver->PostStage(stage,TS->U,solver,mpi,stagetime); CHECKERR(ierr); }
  }

  /* Stage completion */
  for (stage = 0; stage < params->nstages; stage++) {
    ierr = ArrayAXPY(TS->Udot[stage],solver->dt*params->b[stage],solver->u,size*solver->nvars);
    CHECKERR(ierr);
  }

  if (solver->PostStep)  { ierr = solver->PostStep(solver->u,solver,mpi,TS->waqt); CHECKERR(ierr); }

  return(0);
}

