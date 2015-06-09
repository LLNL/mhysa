#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "TransferVecFromPETSc"

int TransferVecFromPETSc(double *u,Vec Y,void *ctxt) 
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = (HyPar*)        context->solver;
  PetscErrorCode  ierr     = 0;
  const double    *Yarr;

  PetscFunctionBegin;

  int *index = (int*) calloc (solver->ndims,sizeof(int));

  ierr = VecGetArrayRead(Y,&Yarr); CHKERRQ(ierr);
  ierr = ArrayCopynD(solver->ndims,(double*)Yarr,u,solver->dim_local,0,
                     solver->ghosts,index,solver->nvars); CHECKERR(ierr);
  ierr = VecRestoreArrayRead(Y,&Yarr); CHKERRQ(ierr);

  free(index);
  PetscFunctionReturn(0);
}

#endif
