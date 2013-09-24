#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <mpivars.h>
#include <petscinterface.h>
#include <hypar.h>

PetscErrorCode PetscPostTimeStep(TS ts)
{
  PETScContext    *context  = NULL;
  HyPar           *solver   = NULL;
  MPIVariables    *mpi      = NULL;
  PetscErrorCode  ierr;
  Vec             Y;

  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscPostTimeStep: Null context!\n");
    return(1);
  }
  solver  = context->solver;
  mpi     = context->mpi;

  int     iter; ierr = TSGetTimeStepNumber(ts,&iter); CHKERRQ(ierr);
  double  dt;   ierr = TSGetTimeStep      (ts,&dt  ); CHKERRQ(ierr);
  double  waqt; ierr = TSGetTime          (ts,&waqt); CHKERRQ(ierr);

  /* get solution */
  ierr = TSGetSolution(ts,&Y); CHKERRQ(ierr);
  ierr = TransferFromPETSc(solver->u,Y,context);

  /* Call any physics-specific post-step function */
  if (solver->PostStep)  { ierr = solver->PostStep(solver->u,solver,mpi,waqt); CHECKERR(ierr); }

  /* Calculate CFL and diffusion number */
  double local_max_cfl  = -1.0, max_cfl  = -1.0;
  double local_max_diff = -1.0, max_diff = -1.0;
  if (solver->ComputeCFL       ) local_max_cfl  = solver->ComputeCFL        (solver,mpi,dt);
  if (solver->ComputeDiffNumber) local_max_diff = solver->ComputeDiffNumber (solver,mpi,dt);
  ierr = MPIMax_double(&max_cfl ,&local_max_cfl ,1); CHECKERR(ierr);
  ierr = MPIMax_double(&max_diff,&local_max_diff,1); CHECKERR(ierr);

  if ((!mpi->rank) && ((iter+1)%solver->screen_op_iter == 0)) {
    printf("Iteration: %6d  "       ,iter+1  );
    printf("Time: %1.3E  "          ,waqt    );
    printf("Max CFL: %1.3E  "       ,max_cfl );
    printf("Max Diff. No.: %1.3E  " ,max_diff);
    printf("\n");
    /* print physics-specific info, if available */
    if (solver->PrintStep) solver->PrintStep(solver,mpi,waqt);
  }

  return(0);
}

#endif
