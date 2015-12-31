/*! @file TransferFromPETSc.c
    @brief Copy from PETSc vector to HyPar array
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "TransferVecFromPETSc"

/*! Copy data from a PETSc vector (used by PETSc time integrators, and with no 
    ghost points) to a HyPar::u array (with ghost points).

    \sa TransferVecToPETSc()
*/
int TransferVecFromPETSc(
                          double  *u,   /*!< HyPar::u array (with ghost points) */
                          Vec     Y,    /*!< PETSc vector */
                          void    *ctxt /*!< Object of type #PETScContext */
                        ) 
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
