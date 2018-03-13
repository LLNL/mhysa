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
  int error = TimeInitialize(solver,mpi,&TS);
  if (error) return(error);

  if (!mpi->rank) printf("Solving in time...\n");
  while(1) {

    /* check for exit conditions */
    if ((TS.t_final >= 0) && (TS.waqt >= TS.t_final)) {
      if (!mpi->rank) printf("Final simulation time (%lf) reached. Exiting time integration.\n",
                             TS.t_final);
      break;
    }
    if ((TS.n_iter >=0) && (TS.iter >= TS.n_iter)) {
      if (!mpi->rank) printf("Maximum number of iterations (%d) reached. Exiting time integration.\n",
                             TS.n_iter);
      break;
    }

    /* Write initial solution to file if this is the first iteration */
    if (!TS.iter) { 
      if (solver->PhysicsOutput) {
        IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
      }
      IERR OutputSolution(solver,mpi); CHECKERR(ierr); 
    }

    /* Call pre-step function */
    IERR TimePreStep  (&TS); CHECKERR(ierr);

    /* Check if this is final step */
    int flag_final_step = 0;
    if (((TS.t_final >= 0) && (TS.waqt + TS.dt >= TS.t_final)) || (TS.iter == (TS.n_iter-1))) {
      flag_final_step = 1;
    }

#ifdef compute_rhs_operators
    /* compute and write (to file) matrix operators representing the right-hand side */
    if (((TS.iter+1)%solver->file_op_iter == 0) || (!TS.iter)) 
      { IERR ComputeRHSOperators(solver,mpi,TS.waqt); CHECKERR(ierr); }
#endif

    /* Step in time */
    IERR TimeStep     (&TS); CHECKERR(ierr);

    /* Call post-step function */
    IERR TimePostStep (&TS); CHECKERR(ierr);

    /* Print information to screen */
    IERR TimePrintStep(&TS, flag_final_step); CHECKERR(ierr);
    tic++;

    /* Write intermediate solution to file */
    if ((TS.iter+1)%solver->file_op_iter == 0) { 
      if (solver->PhysicsOutput) {
        IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
      }
      IERR OutputSolution(solver,mpi); CHECKERR(ierr); 
      tic = 0; 
    }

    TS.iter++;
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
