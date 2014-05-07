#include <stdio.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimePostStep(void *ts)
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
  }

  return(0);
}

