#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <bandedmatrix.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "TransferVecToPETSc"

int TransferVecToPETSc(double *u,Vec Y,void *ctxt) 
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

int TransferMatToPETSc(void *J,Mat A,void *ctxt)
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
}

#endif
