/*! @file PetscPostTimeStep.c
    @brief Post-time-step function
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPostTimeStep"

/*! Function called after every time step */
PetscErrorCode PetscPostTimeStep(TS ts /*!< Time integrator object */)
{
  PETScContext    *context  = NULL;
  HyPar           *solver   = NULL;
  MPIVariables    *mpi      = NULL;
  PetscErrorCode  ierr;
  Vec             Y;

  PetscFunctionBegin;

  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscPostTimeStep: Null context!\n");
    return(1);
  }
  solver  = context->solver;
  mpi     = context->mpi;
  context->tic++;

  int     iter; ierr = TSGetStepNumber(ts,&iter); CHKERRQ(ierr);
  double  dt;   ierr = TSGetTimeStep  (ts,&dt  ); CHKERRQ(ierr);
  double  waqt; ierr = TSGetTime      (ts,&waqt); CHKERRQ(ierr);

  /* get solution */
  ierr = TSGetSolution(ts,&Y); CHKERRQ(ierr);
  ierr = TransferVecFromPETSc(solver->u,Y,context);

  /* Call any physics-specific post-step function */
  if (solver->PostStep)  { ierr = solver->PostStep(solver->u,solver,mpi,waqt,iter); CHECKERR(ierr); }

  /* Calculate CFL and diffusion number */
  double local_max_cfl  = -1.0, max_cfl  = -1.0;
  double local_max_diff = -1.0, max_diff = -1.0;
  if (solver->ComputeCFL       ) local_max_cfl  = solver->ComputeCFL        (solver,mpi,dt,waqt);
  if (solver->ComputeDiffNumber) local_max_diff = solver->ComputeDiffNumber (solver,mpi,dt,waqt);
  ierr = MPIMax_double(&max_cfl ,&local_max_cfl ,1,&mpi->world); CHECKERR(ierr);
  ierr = MPIMax_double(&max_diff,&local_max_diff,1,&mpi->world); CHECKERR(ierr);

  /* Calculate norm of the change in the solution for this time step */
  _ArrayAXPY_(solver->u,-1.0,solver->u0,(solver->npoints_local_wghosts*solver->nvars));
  double sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                                solver->ghosts,solver->index,solver->u0);
  double global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  double norm = sqrt((global_sum/(double)solver->npoints_global));
  
  if (!strcmp(solver->ConservationCheck,"yes")) {
    /* calculate volume integral of the solution at this time step */
    IERR solver->VolumeIntegralFunction(solver->VolumeIntegral,solver->u,solver,mpi); CHECKERR(ierr);
    /* calculate surface integral of the flux at this time step */
    IERR solver->BoundaryIntegralFunction(solver,mpi); CHECKERR(ierr);
    /* calculate the conservation error at this time step       */
    IERR solver->CalculateConservationError(solver,mpi); CHECKERR(ierr);
  }

  if ((!mpi->rank) && (iter%solver->screen_op_iter == 0)) {
    printf("Iteration: %6d  "       ,iter    );
    printf("Time: %1.3E  "          ,waqt    );
    printf("Max CFL: %1.3E  "       ,max_cfl );
    printf("Max Diff. No.: %1.3E  " ,max_diff);
    printf("Norm: %1.4E "           ,norm    );
    /* calculate and print conservation error */
    if (!strcmp(solver->ConservationCheck,"yes")) {
      double error = 0;
      int v;
      for (v=0; v<solver->nvars; v++) 
        error += solver->ConservationError[v] * solver->ConservationError[v];
      error = sqrt(error);
      printf("Conservation loss: %1.4E",error);
    }
    printf("\n");

    /* print physics-specific info, if available */
    if (solver->PrintStep) solver->PrintStep(solver,mpi,waqt);
  }

  /* Write intermediate solution to file */
  if (iter%solver->file_op_iter == 0) { 
    if (solver->PhysicsOutput) {
      IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
    }
    ierr = OutputSolution(solver,mpi); CHECKERR(ierr); context->tic=0; 
  }

  PetscFunctionReturn(0);
}

#endif
