#include <stdio.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

#ifdef compute_rhs_operators
int ComputeRHSOperators(void*,void*,double);
#endif

int Solve(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  /* Define and initialize the time-integration object */
  TimeIntegration TS;
  if (!mpi->rank) printf("Setting up time integration.\n");
  IERR TimeInitialize(solver,mpi,&TS); CHECKERR(ierr);

#ifdef compute_rhs_operators
  /* compute and write (to file) matrix operators representing the right-hand side */
  IERR ComputeRHSOperators(solver,mpi,TS.waqt); CHECKERR(ierr);
#endif
  /* Write an initial solution file */
  IERR OutputSolution(solver,mpi); CHECKERR(ierr);

  if (!mpi->rank) printf("Solving in time (from %d to %d iterations)\n",TS.restart_iter,TS.n_iter);
  for (TS.iter = TS.restart_iter; TS.iter < TS.n_iter; TS.iter++) {
    IERR TimePreStep  (&TS); CHECKERR(ierr);
    IERR TimeStep     (&TS); CHECKERR(ierr);
    IERR TimePostStep (&TS); CHECKERR(ierr);
    IERR TimePrintStep(&TS); CHECKERR(ierr);

#ifdef compute_rhs_operators
    /* compute and write (to file) matrix operators representing the right-hand side */
    if ((TS.iter+1)%solver->file_op_iter == 0) 
      IERR ComputeRHSOperators(solver,mpi,TS.waqt); CHECKERR(ierr);
#endif
    /* Write intermediate solution to file */
    if ((TS.iter+1)%solver->file_op_iter == 0) 
      IERR OutputSolution(solver,mpi); CHECKERR(ierr);
  }

  if (!mpi->rank) printf("Completed time integration (Final time: %f).\n",TS.waqt);
  IERR TimeError(&TS);   CHECKERR(ierr);
  IERR TimeCleanup(&TS); CHECKERR(ierr);
  return(0);
}
