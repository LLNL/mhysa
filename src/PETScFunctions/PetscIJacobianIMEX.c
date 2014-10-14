#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscIJacobianIMEX"

PetscErrorCode PetscIJacobianIMEX(TS ts,PetscReal t,Vec Y,Vec Ydot,PetscReal a,
                                  Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  PetscFunctionBegin;
  context->shift = a;
  context->waqt  = t;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionIMEX"

PetscErrorCode PetscJacobianFunctionIMEX(Mat Jacobian,Vec Y,Vec F)
{
  PETScContext    *context = NULL;
  HyPar           *solver  = NULL;
  MPIVariables    *mpi     = NULL;
  int             ierr     = 0, d;

  PetscFunctionBegin;

  ierr   = MatShellGetContext(Jacobian,&context); CHKERRQ(ierr);
  solver = context->solver;
  mpi    = context->mpi;

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
    /* [J]Y = aY - F(Y) */
    ierr = VecAXPBY(F,context->shift,0,Y);                                CHKERRQ(ierr);

  } else {
    
    double epsilon = 1e-7 / normY;
    /* copy solution from PETSc vector */
    ierr = TransferFromPETSc(u,Y,context);                                CHECKERR(ierr);
    _ArrayAYPX_(uref,epsilon,u,size*solver->nvars);
    /* apply boundary conditions and exchange data over MPI interfaces */
    ierr = solver->ApplyBoundaryConditions(solver,mpi,u   ,NULL,0,t);     CHECKERR(ierr);
    ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                   solver->ghosts,mpi,u);                 CHECKERR(ierr);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
      if (context->flag_hyperbolic_f == _IMPLICIT_) {
        ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->FFunction,solver->Upwind);  
        CHECKERR(ierr);
        _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
        ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF); 
        CHECKERR(ierr);
        _ArrayAXPY_(solver->hyp, 1.0,rhs,size*solver->nvars);
      } 
      if (context->flag_hyperbolic_df == _IMPLICIT_) {
        ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF); 
        CHECKERR(ierr);
        _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
      }
    } else {
      if (context->flag_hyperbolic == _IMPLICIT_) {
        ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->FFunction,solver->Upwind);  
        CHECKERR(ierr);
        _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
      }
    }
    if (context->flag_parabolic == _IMPLICIT_) {
      ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
      _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
    }
    if (context->flag_source == _IMPLICIT_) {
      ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
      _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);
    }

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferToPETSc(rhs,F,context);                  CHECKERR(ierr);
    /* [J]Y = aY - F(Y) */
    ierr = VecAXPBY(F,context->shift,(-1.0/epsilon),Y);     CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#endif
