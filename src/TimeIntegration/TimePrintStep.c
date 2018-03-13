/*! @file TimePrintStep.c
    @brief Print to screen
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Print information to screen (also calls any physics-specific 
  printing function, if defined).
*/
int TimePrintStep(void *ts, /*!< Object of type #TimeIntegration */
                  int flag  /*!< Flag to force printing info */ )
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;
  int             v;

  if ((!mpi->rank) && (((TS->iter+1)%solver->screen_op_iter == 0) || flag)) {
    printf("Iteration: %7d  "   ,TS->iter+1     );
    printf("dt: %1.3E  "        ,TS->dt         );
    printf("Time: %1.3E  "      ,TS->waqt       );
    printf("CFL: %1.3E  "       ,TS->max_cfl    );
    printf("Norm: %1.4E  "      ,TS->norm       );
    printf("wctime: %1.1E (s)  ",TS->iter_wctime);
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

