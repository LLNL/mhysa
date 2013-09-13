#include <stdio.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimePrintStep(void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;

  if ((!mpi->rank) && ((TS->iter+1)%solver->screen_op_iter == 0)) {
    printf("Iteration: %6d  " ,TS->iter+1);
    printf("Time: %E  "       ,TS->waqt);
    printf("Max CFL: %E  "    ,TS->max_cfl);
    printf("Norm: %E  "       ,TS->norm);
    printf("\n");
  }

  return(0);
}

