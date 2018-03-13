/*! @file TimePostStep.c
    @brief Post-time-step function
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Post-time-step function: this function is called at the end of 
  each time step.
  + It updates the current simulation time.
  + It calls functions to print information and to write 
    transient solution to file.
  + It will also call any physics-specific post-time-step function,
    if defined.
*/
int TimePostStep(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;

  /* update current time */
  TS->waqt += TS->dt;

  if ((TS->iter+1)%solver->screen_op_iter == 0) {
    int size = 1, d;
    for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

    /* Calculate norm for this time step */
    _ArrayAXPY_(solver->u,-1.0,TS->u,size*solver->nvars);
    double sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                                  solver->ghosts,solver->index,TS->u);
    double global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    TS->norm = sqrt((global_sum/(double)solver->npoints_global));

    /* write to file */
    if (TS->ResidualFile) fprintf((FILE*)TS->ResidualFile,"%10d\t%E\t%E\n",TS->iter+1,TS->waqt,TS->norm);
  
  }

  if (!strcmp(solver->ConservationCheck,"yes")) {
    /* calculate volume integral of the solution at this time step */
    IERR solver->VolumeIntegralFunction(solver->VolumeIntegral,solver->u,solver,mpi); CHECKERR(ierr);
    /* calculate surface integral of the flux at this time step */
    IERR solver->BoundaryIntegralFunction(solver,mpi); CHECKERR(ierr);
    /* calculate the conservation error at this time step       */
    IERR solver->CalculateConservationError(solver,mpi); CHECKERR(ierr);
  }

  if (solver->PostStep)  { IERR solver->PostStep(solver->u,solver,mpi,TS->waqt,TS->iter); CHECKERR(ierr); }

  gettimeofday(&TS->iter_end_time,NULL);
  long long walltime;
  walltime = (  (TS->iter_end_time.tv_sec * 1000000 + TS->iter_end_time.tv_usec)
              - (TS->iter_start_time.tv_sec * 1000000 + TS->iter_start_time.tv_usec));
  TS->iter_wctime = (double) walltime / 1000000.0;

  return(0);
}

