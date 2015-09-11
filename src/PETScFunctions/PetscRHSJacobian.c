/*! @file PetscRHSJacobian.c
    @brief Contains the functions required for Jacobian computations for implicit time-integration
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscRHSJacobian"
/*! 
    Compute the Jacobian of the right-hand-side \a f(y) for the implicit
    time integration of the ODE \a dy/dt = \a f(y).\n\n
    Matrix-free representation of the Jacobian: This function just saves the  
    time-step size required for the evaluation of the Jacobian. The action of the Jacobian
    are defined through #PetscJacobianFunctionImpl_JFNK (nonlinear problemts) and 
    #PetscJacobianFunctionImpl_Linear (linear problems).\n\n
    The Jacobian matrix is defined as the PETSc type MatShell 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
*/
PetscErrorCode PetscRHSJacobian(
                                 TS ts,        /*!< Time stepping object (see PETSc TS)*/
                                 PetscReal t,  /*!< Current time */
                                 Vec Y,        /*!< Solution vector */
                                 Mat A,        /*!< Jacobian matrix */
                                 Mat B,        /*!< Preconditioning matrix */
                                 void *ctxt    /*!< Application context */
                               )
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;
  _DECLARE_IERR_;

  PetscFunctionBegin;
  solver->count_IJacobian++;
  context->waqt  = t;
  /* Construct preconditioning matrix */
  if (context->flag_use_precon) { IERR PetscComputePreconMatImpl(B,Y,context); CHECKERR(ierr); }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionImpl_JFNK"
/*! 
    This function defines the action of the Jacobian (for implicit 
    time-integration) on a vector for a nonlinear problem. It is evaluated 
    by computing the directional derivative (Jacobian-free Newton-Krylov method)
*/
PetscErrorCode PetscJacobianFunctionImpl_JFNK(
                                              Mat Jacobian, /*!< Jacobian matrix */
                                              Vec Y,        /*!< Input vector */
                                              Vec F         /*!< Output vector (Jacobian times input vector */
                                             )
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
    
    double epsilon =  context->jfnk_eps / normY;
    /* copy solution from PETSc vector */
    ierr = TransferVecFromPETSc(u,Y,context);                             CHECKERR(ierr);
    _ArrayAYPX_(uref,epsilon,u,size*solver->nvars);
    /* apply boundary conditions and exchange data over MPI interfaces */
    ierr = solver->ApplyBoundaryConditions(solver,mpi,u   ,NULL,0,t);     CHECKERR(ierr);
    ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                   solver->ghosts,mpi,u);                 CHECKERR(ierr);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);  CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferVecToPETSc(rhs,F,context);                  CHECKERR(ierr);

  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionImpl_Linear"
/*!
    This function defines the action of the Jacobian (for implicit 
    time-integration) on a vector for a linear problem. 
*/
PetscErrorCode PetscJacobianFunctionImpl_Linear(
                                                Mat Jacobian, /*!< Jacobian matrix */
                                                Vec Y,        /*!< Input vector */
                                                Vec F         /*!< Output vector (Jacobian times input vector */
                                               )
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
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);  CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferVecToPETSc(rhs,F,context);     CHECKERR(ierr);
  }

  PetscFunctionReturn(0);
}

#endif
