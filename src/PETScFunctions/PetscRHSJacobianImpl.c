#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

/* Unpreconditioned Jacobian-free Newton Krylov */
#undef __FUNCT__
#define __FUNCT__ "PetscRHSJacobianImpl_JFNK_NoPre"
PetscErrorCode PetscRHSJacobianImpl_JFNK_NoPre(TS ts,PetscReal t,Vec Y,Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;

  PetscFunctionBegin;
  solver->count_RHSJacobian++;
  context->waqt  = t;
  PetscFunctionReturn(0);
}

/* Preconditioned Jacobian-free Newton Krylov */
#undef __FUNCT__
#define __FUNCT__ "PetscRHSJacobianImpl_JFNK_Pre"
PetscErrorCode PetscRHSJacobianImpl_JFNK_Pre(TS ts,PetscReal t,Vec Y,Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;
  MPIVariables *mpi     = context->mpi;
  _DECLARE_IERR_;
  PetscFunctionBegin;
  solver->count_RHSJacobian++;

  double *u = solver->u;
  void   *J = solver->Jac;

  /* copy solution from PETSc vector */
  IERR TransferVecFromPETSc(u,Y,context);                         CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);    CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,u);             CHECKERR(ierr);

  /* evaluate the Jacobian function */
  IERR solver->PFunction(J,u,solver,mpi,t,0);                     CHECKERR(ierr);
  /* transfer the Jacobian to PETSc Mat */
  IERR TransferMatToPETSc(J,B,context);                           CHECKERR(ierr);

  context->waqt  = t;
  PetscFunctionReturn(0);
}

/* Specified Jacobian and preconditioning matrices */
#undef __FUNCT__
#define __FUNCT__ "PetscRHSJacobianImpl_Jac_Pre"
PetscErrorCode PetscRHSJacobianImpl_Jac_Pre(TS ts,PetscReal t,Vec Y,Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;
  MPIVariables *mpi     = context->mpi;
  _DECLARE_IERR_;
  PetscFunctionBegin;
  solver->count_RHSJacobian++;

  double *u = solver->u;
  void   *J = solver->Jac;

  /* copy solution from PETSc vector */
  IERR TransferVecFromPETSc(u,Y,context);                         CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);    CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,u);             CHECKERR(ierr);

  /* evaluate the Jacobian function */
  IERR solver->JFunction(J,u,solver,mpi,t,0);                     CHECKERR(ierr);
  /* transfer the Jacobian to PETSc Mat */
  IERR TransferMatToPETSc(J,A,context);                           CHECKERR(ierr);

  /* evaluate the preconditioning function */
  IERR solver->PFunction(J,u,solver,mpi,t,0);                     CHECKERR(ierr);
  /* transfer the Jacobian to PETSc Mat */
  IERR TransferMatToPETSc(J,B,context);                           CHECKERR(ierr);
  
  /* Done */
  PetscFunctionReturn(0);
}

/* Jacobian-free Newton Krylov with specified Jacobian matrix as preconditioner */
#undef __FUNCT__
#define __FUNCT__ "PetscRHSJacobianImpl_JFNK_JacIsPre"
PetscErrorCode PetscRHSJacobianImpl_JFNK_JacIsPre(TS ts,PetscReal t,Vec Y,Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;
  MPIVariables *mpi     = context->mpi;
  _DECLARE_IERR_;
  PetscFunctionBegin;
  solver->count_RHSJacobian++;

  double *u = solver->u;
  void   *J = solver->Jac;

  /* copy solution from PETSc vector */
  IERR TransferVecFromPETSc(u,Y,context);                         CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);    CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,u);             CHECKERR(ierr);

  /* evaluate the Jacobian function */
  IERR solver->JFunction(J,u,solver,mpi,t,0);                     CHECKERR(ierr);

  /* transfer the Jacobian to PETSc Mat */
  IERR TransferMatToPETSc(J,B,context);                           CHECKERR(ierr);

  context->waqt  = t;
  /* Done */
  PetscFunctionReturn(0);
}

/* Specified Jacobian matrix, used as preconditioner as well */
#undef __FUNCT__
#define __FUNCT__ "PetscRHSJacobianImpl_Jac_NoPre"
PetscErrorCode PetscRHSJacobianImpl_Jac_NoPre(TS ts,PetscReal t,Vec Y,Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;
  MPIVariables *mpi     = context->mpi;
  _DECLARE_IERR_;
  PetscFunctionBegin;
  solver->count_RHSJacobian++;

  /* check that A & B are the same matrices */
  if (A != B) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in PetscRHSJacobianImpl_Jac_NoPre(): ");
      fprintf(stderr,"A and B need to be the same Mat object.\n");
    }
    PetscFunctionReturn(1);
  }

  double *u = solver->u;
  void   *J = solver->Jac;

  /* copy solution from PETSc vector */
  IERR TransferVecFromPETSc(u,Y,context);                         CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);    CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,u);             CHECKERR(ierr);

  /* evaluate the Jacobian function */
  IERR solver->JFunction(J,u,solver,mpi,t,0);                     CHECKERR(ierr);

  /* transfer the Jacobian to PETSc Mat */
  IERR TransferMatToPETSc(J,A,context);                           CHECKERR(ierr);

  /* Done */
  PetscFunctionReturn(0);
}

/* Directional derivative evaluation of Jacobian times a vector for the
   Jacobian-free Newton-Krylov approach */
#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionImpl_JFNK"
PetscErrorCode PetscJacobianFunctionImpl_JFNK(Mat Jacobian,Vec Y,Vec F)
{
  PETScContext    *context = NULL;
  HyPar           *solver  = NULL;
  MPIVariables    *mpi     = NULL;
  int             ierr     = 0, d;

  PetscFunctionBegin;

  ierr   = MatShellGetContext(Jacobian,&context); CHKERRQ(ierr);
  solver = context->solver;
  mpi    = context->mpi;
  solver->count_IJacFunction++;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u       = solver->u;
  double *uref    = solver->uref;
  double *rhsref  = solver->rhsref;
  double *rhs     = solver->rhs;

  double t = context->waqt; /* current stage/step time */

  double normY;
  ierr = VecNorm(Y,NORM_2,&normY);                                        CHKERRQ(ierr);

  if (normY < 1e-16) {

    /* F = 0 */
    ierr = VecZeroEntries(F);                                             CHKERRQ(ierr);

  } else {
    
    double epsilon = (context->flag_is_linear ? 1.0 : 1e-7 / normY );
    /* copy solution from PETSc vector */
    ierr = TransferVecFromPETSc(u,Y,context);                             CHECKERR(ierr);
    _ArrayAYPX_(uref,epsilon,u,size*solver->nvars);
    /* apply boundary conditions and exchange data over MPI interfaces */
    ierr = solver->ApplyBoundaryConditions(solver,mpi,u   ,NULL,0,t);     CHECKERR(ierr);
    ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                   solver->ghosts,mpi,u);                 CHECKERR(ierr);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);  CHECKERR(ierr);
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    _ArrayAXPY_(solver->hyp,   -1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->par,    1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferVecToPETSc(rhs,F,context);                             CHECKERR(ierr);
  }

  PetscFunctionReturn(0);
}

/* Function to evaluate Jacobian times a vector for a linear case */
#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionImpl_Linear"
PetscErrorCode PetscJacobianFunctionImpl_Linear(Mat Jacobian,Vec Y,Vec F)
{
  PETScContext    *context = NULL;
  HyPar           *solver  = NULL;
  MPIVariables    *mpi     = NULL;
  int             ierr     = 0, d;

  PetscFunctionBegin;

  ierr   = MatShellGetContext(Jacobian,&context); CHKERRQ(ierr);
  solver = context->solver;
  mpi    = context->mpi;
  solver->count_IJacFunction++;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u       = solver->u;
  double *uref    = solver->uref;
  double *rhsref  = solver->rhsref;
  double *rhs     = solver->rhs;

  double t = context->waqt; /* current stage/step time */

  double normY;
  ierr = VecNorm(Y,NORM_2,&normY);                                        CHKERRQ(ierr);

  if (normY < 1e-16) {

    /* F = 0 */
    ierr = VecZeroEntries(F);                                             CHKERRQ(ierr);

  } else {
    
    /* copy solution from PETSc vector */
    ierr = TransferVecFromPETSc(u,Y,context);                             CHECKERR(ierr);
    _ArrayAYPX_(uref,1.0,u,size*solver->nvars);
    /* apply boundary conditions and exchange data over MPI interfaces */
    ierr = solver->ApplyBoundaryConditions(solver,mpi,u   ,NULL,0,t);     CHECKERR(ierr);
    ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                   solver->ghosts,mpi,u);                 CHECKERR(ierr);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);  CHECKERR(ierr);
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    _ArrayAXPY_(solver->hyp   ,-1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->par   , 1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferVecToPETSc(rhs,F,context);                             CHECKERR(ierr);
  }

  PetscFunctionReturn(0);
}

#endif
