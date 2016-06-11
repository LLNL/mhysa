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
  int npoints = context->npoints, 
      *points = context->points,
      nvars   = solver->nvars,
      ndims   = solver->ndims,
      nv      = ndims + 1, n;

  ierr = VecGetArrayRead(Y,&Yarr); CHKERRQ(ierr);
  for (n = 0; n < npoints; n++) {
    int p = (points+n*nv)[ndims];
    _ArrayCopy1D_((Yarr+n*nvars),(u+p*nvars),nvars);
  }
  ierr = VecRestoreArrayRead(Y,&Yarr); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif
