/*! @file PetscRHSFunctionExpl.c
    @brief Function to compute the right-hand-side for explicit time integration
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
#define __FUNCT__ "PetscRHSFunctionExpl"

/*!
  Compute the right-hand-side (RHS) for the explicit time integration of the 
  governing equations: The spatially discretized ODE can be expressed as
  \f{equation}{
    \frac {d{\bf U}} {dt} = {\bf F}\left({\bf U}\right).
  \f}
  This function computes \f${\bf F}\left({\bf U}\right)\f$, given \f${\bf U}\f$.

  \sa TSSetRHSFunction (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSFunction.html)

  Note:
  + \f${\bf U}\f$ is \a Y in the code.
  + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html for documentation on
    PETSc's time integrators.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscRHSFunctionExpl(
                                      TS        ts,   /*!< Time integration object */
                                      PetscReal t,    /*!< Current simulation time */
                                      Vec       Y,    /*!< State vector (input) */
                                      Vec       F,    /*!< The computed right-hand-side vector */
                                      void      *ctxt /*!< Object of type #PETScContext */
                                   )
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = (HyPar*)        context->solver;
  MPIVariables    *mpi     = (MPIVariables*) context->mpi;
  int             d;
  _DECLARE_IERR_;

  PetscFunctionBegin;
  solver->count_RHSFunction++;
  
  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u   = solver->u;
  double *rhs = solver->rhs;

  /* copy solution from PETSc vector */
  IERR TransferVecFromPETSc(u,Y,context);                                         CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);                      CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,u);                             CHECKERR(ierr);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->FFunction,solver->Upwind);
                                                                                  CHECKERR(ierr);
  IERR solver->ParabolicFunction (solver->par,u,solver,mpi,t);                    CHECKERR(ierr);
  IERR solver->SourceFunction    (solver->source,u,solver,mpi,t);                 CHECKERR(ierr);

  _ArraySetValue_(rhs,size*solver->nvars,0.0);
  _ArrayAXPY_(solver->hyp   ,-1.0,rhs,size*solver->nvars);
  _ArrayAXPY_(solver->par   , 1.0,rhs,size*solver->nvars);
  _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

  /* Transfer RHS to PETSc vector */
  IERR TransferVecToPETSc(rhs,F,context);                                         CHECKERR(ierr);

  PetscFunctionReturn(0);
}

#endif
