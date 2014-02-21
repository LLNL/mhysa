#include <stdio.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int Solve(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  /* Define and initialize the time-integration object */
  TimeIntegration TS;
  if (!mpi->rank) printf("Setting up time integration.\n");
  IERR TimeInitialize(solver,mpi,&TS); CHECKERR(ierr);

  if (!mpi->rank) printf("Solving in time (%d iterations)\n",TS.n_iter);
  for (TS.iter = 0; TS.iter < TS.n_iter; TS.iter++) {
    IERR TimePreStep  (&TS); CHECKERR(ierr);
    IERR TimeStep     (&TS); CHECKERR(ierr);
    IERR TimePostStep (&TS); CHECKERR(ierr);
    IERR TimePrintStep(&TS); CHECKERR(ierr);

    /* Write intermediate solution to file */
    if ((TS.iter+1)%solver->file_op_iter == 0) 
      IERR OutputSolution(solver,mpi); CHECKERR(ierr);
  }

  if (!mpi->rank) printf("Completed time integration (Final time: %f).\n",TS.waqt);
  IERR TimeCleanup(&TS); CHECKERR(ierr);
  return(0);
}
