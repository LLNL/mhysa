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
  _DECLARE_IERR_;

  /* copy current solution for norm computation later */
  if ((TS->iter+1)%solver->screen_op_iter == 0) {
    int size = 1,d;
    for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d] + 2*solver->ghosts);
    _ArrayCopy1D_(solver->u,TS->u,size*solver->nvars);

    /* compute max CFL and diffusion number over the domain */
    double local_max_cfl  = -1.0;
    double local_max_diff = -1.0;
    if (solver->ComputeCFL       ) local_max_cfl  = solver->ComputeCFL        (solver,mpi,TS->dt,TS->waqt);
    if (solver->ComputeDiffNumber) local_max_diff = solver->ComputeDiffNumber (solver,mpi,TS->dt,TS->waqt);
    IERR MPIMax_double(&TS->max_cfl ,&local_max_cfl ,1,&mpi->world); CHECKERR(ierr);
    IERR MPIMax_double(&TS->max_diff,&local_max_diff,1,&mpi->world); CHECKERR(ierr);
  }

  /* set the step boundary flux integral value to zero */
  _ArraySetValue_(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);

  return(0);
}

