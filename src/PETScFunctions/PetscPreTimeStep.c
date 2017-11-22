/*! @file PetscPreTimeStep.c
    @brief Pre-time-step function
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <petscinterface.h>
#include <hypar.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPreTimeStep"

#ifdef compute_rhs_operators
int PetscComputeRHSOperators(TS,double,void*);
#endif

/*! Function called before a time step */
PetscErrorCode PetscPreTimeStep(TS ts /*!< Time integration object */)
{
  PETScContext    *context  = NULL;
  HyPar           *solver   = NULL;
  MPIVariables    *mpi      = NULL;
  PetscErrorCode  ierr      = 0;
  Vec             Y;
  TSType          time_scheme;
  double          waqt;
  int             iter;

  PetscFunctionBegin;

  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscPreTimeStep: Null context!\n");
    return(1);
  }
  solver  = context->solver;
  mpi     = context->mpi;

  /* get solution */
  ierr = TSGetSolution(ts,&Y);                        CHKERRQ(ierr);
  ierr = TransferVecFromPETSc(solver->u,Y,context);   CHECKERR(ierr);
  ierr = TSGetTime(ts,&waqt);                         CHKERRQ(ierr);
  ierr = TSGetStepNumber(ts,&iter);                   CHKERRQ(ierr);

  /* save a copy of the solution to compute norm at end of time step */
  _ArrayCopy1D_(solver->u,solver->u0,(solver->npoints_local_wghosts*solver->nvars));

  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,solver->u,NULL,waqt); CHECKERR(ierr);
  IERR solver->ApplyIBConditions(solver,mpi,solver->u,waqt); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,solver->u);             CHECKERR(ierr);

  /* Call any physics-specific pre-step function */
  if (solver->PreStep) { ierr = solver->PreStep(solver->u,solver,mpi,waqt); CHECKERR(ierr); }

  /* If using a non-linear scheme with ARKIMEX methods, 
     compute the non-linear finite-difference operator */
  ierr = TSGetType(ts,&time_scheme);        CHKERRQ(ierr);
  if (!strcmp(time_scheme,TSARKIMEX)) {
    ierr = solver->NonlinearInterp(solver->u,solver,mpi,waqt,solver->FFunction); 
    CHECKERR(ierr);
  }
  
  /* set the step boundary flux integral value to zero */
  _ArraySetValue_(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);

#ifdef compute_rhs_operators
  if ((!iter) || ((iter+1)%solver->file_op_iter == 0)) 
    { ierr = PetscComputeRHSOperators(ts,waqt,context); CHECKERR(ierr); }
#endif
  /* Write initial solution file if this is the first iteration */
  if (!iter) { 
    if (solver->PhysicsOutput) {
      IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
    }
    ierr = OutputSolution(solver,mpi); CHECKERR(ierr); 
  }

  IERR TransferVecToPETSc(solver->u,Y,context);
  PetscFunctionReturn(0);
}

#endif
