/*! @file Solve.c
    @author Debojyoti Ghosh
    @brief  Solve the governing equations in time
*/

#include <stdio.h>
#include <mpivars.h>
#include <io.h>
#include <hypar.h>
#include <timeintegration.h>

#ifdef compute_rhs_operators
int ComputeRHSOperators(void*,void*,double);
#endif

/*! This function integrates the semi-discrete ODE (obtained from discretizing the 
    PDE in space) using natively implemented time integration methods. It initializes 
    the time integration object, iterates the simulation for the required number of 
    time steps, and calculates the errors. After the specified number of iterations, 
    it writes out some information to the screen and the solution to a file.
*/
int Solve(
            void *s,  /*!< Solver object of type #HyPar */
            void *m   /*!< MPI object of type #MPIVariables */
         )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           tic     = 0;
  _DECLARE_IERR_;

  /* write out iblank to file for visualization */
  if (solver->flag_ib) {
    WriteArray(solver->ndims,1,solver->dim_global,solver->dim_local,
               solver->ghosts,solver->x,solver->iblank,solver,mpi,"iblank");
  }

  /* Define and initialize the time-integration object */
  TimeIntegration TS;
  if (!mpi->rank) printf("Setting up time integration.\n");
  IERR TimeInitialize(solver,mpi,&TS); CHECKERR(ierr);

  if (!mpi->rank) printf("Solving in time (from %d to %d iterations)\n",TS.restart_iter,TS.n_iter);
  for (TS.iter = TS.restart_iter; TS.iter < TS.n_iter; TS.iter++) {
    /* Call pre-step function */
    IERR TimePreStep  (&TS); CHECKERR(ierr);
#ifdef compute_rhs_operators
    /* compute and write (to file) matrix operators representing the right-hand side */
    if (((TS.iter+1)%solver->file_op_iter == 0) || (!TS.iter)) 
      { IERR ComputeRHSOperators(solver,mpi,TS.waqt); CHECKERR(ierr); }
#endif
    /* Write initial solution to file if this is the first iteration */
    if (!TS.iter) { 
      if (solver->PhysicsOutput) {
        IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
      }
      IERR OutputSolution(solver,mpi); CHECKERR(ierr); 
    }
    /* Step in time */
    IERR TimeStep     (&TS); CHECKERR(ierr);
    /* Call post-step function */
    IERR TimePostStep (&TS); CHECKERR(ierr);
    /* Print information to screen */
    IERR TimePrintStep(&TS); CHECKERR(ierr);
    tic++;

    /* Write intermediate solution to file */
    if ((TS.iter+1)%solver->file_op_iter == 0) { 
      if (solver->PhysicsOutput) {
        IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
      }
      IERR OutputSolution(solver,mpi); CHECKERR(ierr); 
      tic = 0; 
    }
  }

  /* write a final solution file, if last iteration did not write one */
  if (tic || (!TS.n_iter)) { 
    if (solver->PhysicsOutput) {
      IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
    }
    IERR OutputSolution(solver,mpi); CHECKERR(ierr); 
  }

  if (!mpi->rank) printf("Completed time integration (Final time: %f).\n",TS.waqt);

  /* calculate error if exact solution has been provided */
  IERR CalculateError(solver,mpi); CHECKERR(ierr);
  IERR TimeCleanup(&TS); CHECKERR(ierr);
  return(0);
}
