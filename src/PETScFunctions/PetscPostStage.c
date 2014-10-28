#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <mpivars.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPostStage"

PetscErrorCode PetscPostStage(TS ts,PetscReal stagetime,PetscInt stageindex,Vec *Y)
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

  /* If using a non-linear scheme with ARKIMEX methods, 
     compute the non-linear finite-difference operator */
  ierr = TSGetType(ts,&time_scheme); CHKERRQ(ierr);
  if (!strcmp(time_scheme,TSARKIMEX)) {
    ierr = solver->NonlinearInterp(solver->u,solver,mpi,(double)stagetime,
                                   solver->FFunction); CHECKERR(ierr);
  }

  PetscFunctionReturn(0);
}

#endif
