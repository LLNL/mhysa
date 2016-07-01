/*! @file PetscPostStage.c
    @author Debojyoti Ghosh
    @brief Post-time-integration-stage function.
*/
#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <mpivars.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPostStage"

/*! Function called after every stage in a multi-stage time-integration method */
PetscErrorCode PetscPostStage(
                                TS        ts,         /*!< Time integrator of PETSc type TS */
                                PetscReal stagetime,  /*!< Current stage time */
                                PetscInt  stageindex, /*!< Stage */
                                Vec       *Y          /*!< Stage solutions (all stages) - 
                                                           be careful what you access */
                             )
{
  PETScContext    *context  = NULL;
  HyPar           *solver   = NULL;
  MPIVariables    *mpi      = NULL;
  TSType          time_scheme;
  PetscErrorCode  ierr      = 0;

  PetscFunctionBegin;

  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscPreTimeStep: Null context!\n");
    return(1);
  }
  solver  = context->solver;
  mpi     = context->mpi;

  /* get solution */
  ierr = TransferVecFromPETSc(solver->u,Y[stageindex],context); CHECKERR(ierr);

  /* apply immersed boundaries */
  IERR solver->ApplyIBConditions(solver,mpi,solver->u,stagetime); CHECKERR(ierr);

  /* If using a non-linear scheme with ARKIMEX methods, 
     compute the non-linear finite-difference operator */
  ierr = TSGetType(ts,&time_scheme); CHKERRQ(ierr);
  if (!strcmp(time_scheme,TSARKIMEX)) {
    ierr = solver->NonlinearInterp(solver->u,solver,mpi,(double)stagetime,
                                   solver->FFunction); CHECKERR(ierr);
  }

  /* Call any physics-specific post-stage function, if available */
  if (solver->PostStage) {
    ierr = solver->PostStage(solver->u,solver,mpi,stagetime); CHECKERR(ierr);
  }

  IERR TransferVecToPETSc(solver->u,Y[stageindex],context);
  PetscFunctionReturn(0);
}

#endif
