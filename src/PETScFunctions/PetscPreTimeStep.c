#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <mpivars.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPreTimeStep"

PetscErrorCode PetscPreTimeStep(TS ts)
{
  PETScContext    *context  = NULL;
  HyPar           *solver   = NULL;
  MPIVariables    *mpi      = NULL;
  PetscErrorCode  ierr      = 0;
  Vec             Y;
  TSType          time_scheme;
  double          waqt;

  PetscFunctionBegin;

  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscPreTimeStep: Null context!\n");
    return(1);
  }
  solver  = context->solver;
  mpi     = context->mpi;

  /* get solution */
  ierr = TSGetSolution(ts,&Y);                    CHKERRQ(ierr);
  ierr = TransferFromPETSc(solver->u,Y,context);  CHECKERR(ierr);
  ierr = TSGetTime(ts,&waqt);                     CHKERRQ(ierr);

  /* Call any physics-specific pre-step function */
  if (solver->PreStep)  { 
    ierr = solver->PreStep(solver->u,solver,mpi,waqt); 
    CHECKERR(ierr); 
  }

  /* If using a non-linear scheme with ARKIMEX methods, 
     compute the non-linear finite-difference operator */
  ierr = TSGetType(ts,&time_scheme);        CHKERRQ(ierr);
  if (!strcmp(time_scheme,TSARKIMEX)) {
    ierr = solver->NonlinearInterp(solver->u,solver,mpi,waqt,solver->FFunction); 
    CHECKERR(ierr);
  }

  PetscFunctionReturn(0);
}

#endif
