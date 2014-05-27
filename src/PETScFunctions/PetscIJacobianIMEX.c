#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

PetscErrorCode PetscIJacobianIMEX(TS ts,PetscReal t,Vec Y,Vec Ydot,PetscReal a,
                                  Mat A,Mat B,void *ctxt)
{
  PETScContext *context = (PETScContext*) ctxt;
  context->shift = a;
  context->waqt  = t;
  return(0);
}

PetscErrorCode PetscJacobianFunctionIMEX(Mat Jacobian,Vec Y,Vec F)
{
  PETScContext    *context = NULL;
  HyPar           *solver  = NULL;
  MPIVariables    *mpi     = NULL;
  int             ierr     = 0, d;

  ierr   = MatShellGetContext(Jacobian,&context); CHKERRQ(ierr);
  solver = context->solver;
  mpi    = context->mpi;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u    = solver->u;
  double *uref = solver->uref;
  double *rhs  = (double*) calloc (size*solver->nvars,sizeof(double));

  double t = context->waqt; /* current stage/step time */

  /* copy solution from PETSc vector */
  ierr = TransferFromPETSc(u,Y,context);                              CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,u,uref,1,t);      CHECKERR(ierr);
  ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,u);               CHECKERR(ierr);

  /* initialize right-hand side to zero */
  _ArraySetValue_(rhs,size*solver->nvars,0.0);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  if (solver->HyperbolicFunction && (context->flag_hyperbolic == _IMPLICIT_)) {
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t);    CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
  }
  if (solver->HyperbolicFunction1 && (context->flag_hyperbolic1 == _IMPLICIT_)) {
    ierr = solver->HyperbolicFunction1(solver->hyp,u,solver,mpi,t);   CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
  }
  if (solver->HyperbolicFunction2 && (context->flag_hyperbolic2 == _IMPLICIT_)) {
    ierr = solver->HyperbolicFunction2(solver->hyp,u,solver,mpi,t);   CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
  }
  if (solver->ParabolicFunction && (context->flag_parabolic == _IMPLICIT_)) {
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);    CHECKERR(ierr);
    _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
  }
  if (solver->SourceFunction && (context->flag_source == _IMPLICIT_)) {
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t); CHECKERR(ierr);
    _ArrayAXPY_(solver->source , 1.0,rhs,size*solver->nvars);
  }

  /* Transfer RHS to PETSc vector */
  ierr = TransferToPETSc(rhs,F,context);                              CHECKERR(ierr);

  /* [J]Y = aY - F(Y) */
  ierr = VecAXPBY(F,context->shift,-1.0,Y); CHKERRQ(ierr); CHKERRQ(ierr);

  return(0);
}

#endif
