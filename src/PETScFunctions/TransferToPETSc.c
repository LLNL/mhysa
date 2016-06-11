/*! @file TransferToPETSc.c
    @brief Copy from a HyPar array to a PETSc vector.
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <bandedmatrix.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "TransferVecToPETSc"

/*! Copy data to a PETSc vector (used by PETSc time integrators, and with no 
    ghost points) from a HyPar::u array (with ghost points).

    \sa TransferVecFromPETSc()
*/
int TransferVecToPETSc(
                        double  *u,   /*!< HyPar::u array (with ghost points) */
                        Vec     Y,    /*!< PETSc vector */
                        void    *ctxt /*!< Object of type #PETScContext */
                      )
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = (HyPar*)        context->solver;
  PetscErrorCode  ierr     = 0;
  double          *Yarr;

  PetscFunctionBegin;
  int npoints = context->npoints, 
      *points = context->points,
      nvars   = solver->nvars,
      ndims   = solver->ndims,
      nv      = ndims + 1, n;

  ierr = VecGetArray(Y,&Yarr); CHKERRQ(ierr);
  for (n = 0; n < npoints; n++) {
    int p = (points+n*nv)[ndims];
    _ArrayCopy1D_((u+p*nvars),(Yarr+n*nvars),nvars);
  }
  ierr = VecRestoreArray(Y,&Yarr); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*!
  Copy a matrix of type #BandedMatrix to a PETSc matrix.
*/
int TransferMatToPETSc(
                        void *J,    /*!< Matrix of type #BandedMatrix */
                        Mat   A,    /*!< PETSc matrix */
                        void *ctxt  /*!< Object of type #PETScContext */
                      )
{
  BandedMatrix    *M = (BandedMatrix*) J;
  PetscErrorCode  ierr     = 0;
  int             i, n, p, q, bs = M->BlockSize, nbands = M->nbands, bs2 = bs*bs;

  for (i=0; i<M->nrows_local; i++) {
    int     colind[nbands];
    double  val[bs][bs*nbands];
    for (n=0; n<nbands; n++) {
      colind[n] = M->ncol[nbands*i+n];
      for (p=0; p<bs; p++) {
        for (q = 0; q<bs; q++) {
          val[p][n*bs+q] = M->data[i*nbands*bs2+n*bs2+p*bs+q];
        }
      }
    }
    ierr = MatSetValuesBlocked(A,1,&M->nrow[i],M->nbands,&colind[0],&val[0][0],INSERT_VALUES); 
                                                                                  CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);                                  CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);                                  CHKERRQ(ierr);

  return(0);
}

#endif
