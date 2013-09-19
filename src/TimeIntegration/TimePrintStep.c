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
    printf("Iteration: %6d  "       ,TS->iter+1  );
    printf("Time: %1.3E  "          ,TS->waqt    );
    printf("Max CFL: %1.3E  "       ,TS->max_cfl );
    printf("Max Diff. No.: %1.3E  " ,TS->max_diff);
    printf("Norm: %1.4E  "          ,TS->norm    );
    printf("\n");
  }

  return(0);
}

