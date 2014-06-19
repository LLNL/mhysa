#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "TransferToPETSc"

int TransferToPETSc(double *u,Vec Y,void *ctxt) 
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = (HyPar*)        context->solver;
  PetscErrorCode  ierr     = 0;
  double          *Yarr;

  PetscFunctionBegin;

  int *index = (int*) calloc (solver->ndims,sizeof(int));
  ierr = VecGetArray(Y,&Yarr); CHKERRQ(ierr);
  ierr = ArrayCopynD(solver->ndims,u,Yarr,solver->dim_local,solver->ghosts,
                     0,index,solver->nvars); CHECKERR(ierr);
  ierr = VecRestoreArray(Y,&Yarr); CHKERRQ(ierr);

  free(index);
  PetscFunctionReturn(0);
}

#endif
