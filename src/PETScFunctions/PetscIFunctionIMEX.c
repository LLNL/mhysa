#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

PetscErrorCode PetscIFunctionIMEX(TS ts, PetscReal t, Vec Y, Vec Ydot, Vec F, void *ctxt)
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = (HyPar*)        context->solver;
  MPIVariables    *mpi     = (MPIVariables*) context->mpi;
  int             ierr     = 0, d;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u   = solver->u;
  double *rhs = solver->rhs;

  /* copy solution from PETSc vector */
  ierr = TransferFromPETSc(u,Y,context);                            CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);    CHECKERR(ierr);
  ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,u);             CHECKERR(ierr);

  /* initialize right-hand side to zero */
  _ArraySetValue_(rhs,size*solver->nvars,0.0);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (context->flag_hyperbolic_f == _IMPLICIT_) {
      ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->FFunction);  CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
      ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->dFFunction); CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp, 1.0,rhs,size*solver->nvars);
    } 
    if (context->flag_hyperbolic_df == _IMPLICIT_) {
      ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->dFFunction); CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    }
  } else {
    if (context->flag_hyperbolic == _IMPLICIT_) {
      ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->FFunction);  CHECKERR(ierr);
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

  /* save a copy of the solution and RHS for use in IJacobian */
  _ArrayCopy1D_(u  ,solver->uref  ,(size*solver->nvars));
  _ArrayCopy1D_(rhs,solver->rhsref,(size*solver->nvars));

  /* Transfer RHS to PETSc vector */
  ierr = TransferToPETSc(rhs,F,context);                                  CHECKERR(ierr);

  /* LHS = Ydot - F(u) */
  ierr = VecAYPX(F,-1.0,Ydot); CHKERRQ(ierr);

  return(0);
}

#endif
