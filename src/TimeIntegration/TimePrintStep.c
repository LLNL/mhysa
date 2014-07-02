#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimePrintStep(void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;
  int             v;

  if ((!mpi->rank) && ((TS->iter+1)%solver->screen_op_iter == 0)) {
    printf("Iteration: %7d  "       ,TS->iter+1  );
    printf("Time: %1.3E  "          ,TS->waqt    );
    printf("Max CFL: %1.3E  "       ,TS->max_cfl );
    printf("Max Diff. No.: %1.3E  " ,TS->max_diff);
    printf("Norm: %1.4E  "          ,TS->norm    );
    /* calculate and print conservation error */
    if (!strcmp(solver->ConservationCheck,"yes")) {
      double error = 0;
      for (v=0; v<solver->nvars; v++) 
        error += solver->ConservationError[v] * solver->ConservationError[v];
      error = sqrt(error);
      printf("Conservation loss: %1.4E",error);
    }
    printf("\n");
    /* print physics-specific info, if available */
    if (solver->PrintStep) solver->PrintStep(solver,mpi,TS->waqt);
  }

  return(0);
}

