#include <stdio.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int Solve(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0;

  /* Define and initialize the time-integration object */
  TimeIntegration TS;
  ierr = TimeInitialize(solver,mpi,&TS); CHECKERR(ierr);

  if (!mpi->rank) printf("Solving in time (%d iterations)\n",TS.n_iter);
  for (TS.iter = 0; TS.iter < TS.n_iter; TS.iter++) {
    ierr = TimePreStep  (&TS); CHECKERR(ierr);
    ierr = TimeStep     (&TS); CHECKERR(ierr);
    ierr = TimePostStep (&TS); CHECKERR(ierr);
    ierr = TimePrintStep(&TS); CHECKERR(ierr);

    /* Write intermediate solution to file */
    if ((TS.iter+1)%solver->file_op_iter == 0) 
      ierr = OutputSolution(solver,mpi); CHECKERR(ierr);
  }

  if (!mpi->rank) printf("Completed time integration (Final time: %f).\n",TS.waqt);
  ierr = TimeCleanup(&TS); CHECKERR(ierr);
  return(0);
}
