#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <petscinterface.h>
#include <hypar.h>

int TransferFromPETSc(double *u,Vec Y,void *ctxt) 
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = (HyPar*)        context->solver;
  PetscErrorCode  ierr     = 0;
  double          *Yarr;

  int *index = (int*) calloc (solver->ndims,sizeof(int));

  ierr = VecGetArray(Y,&Yarr); CHKERRQ(ierr);
  ierr = ArrayCopynD(solver->ndims,Yarr,u,solver->dim_local,0,
                     solver->ghosts,index,solver->nvars); CHECKERR(ierr);
  ierr = VecRestoreArray(Y,&Yarr); CHKERRQ(ierr);

  free(index);
  return(0);
}

#endif
