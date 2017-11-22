/*! @file PetscRHSFunctionIMEX.c
    @brief Compute the explicitly-treated part of the right-hand side for IMEX time integration
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
#define __FUNCT__ "PetscRHSFunctionIMEX"

/*!
  Compute the explicitly-treated right-hand-side (RHS) for the implicit-explicit (IMEX) time integration
  of the governing equations: The ODE, obtained after discretizing the governing PDE in space,
  is expressed as follows (for the purpose of IMEX time integration):
  \f{eqnarray}{
    \frac {d{\bf U}}{dt} &=& {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right), \\
    \Rightarrow \dot{\bf U} - {\bf G}\left({\bf U}\right) &=& {\bf F}\left({\bf U}\right), 
  \f}
  where \f${\bf F}\f$ is non-stiff and integrated in time explicitly, and \f${\bf G}\f$
  is stiff and integrated in time implicitly, and \f${\bf U}\f$ represents the entire
  solution vector (state vector).

  Note that \f${\bf F}\left({\bf U}\right)\f$ represents all terms that the user has indicated to be
  integrated in time explicitly (#PETScContext::flag_hyperbolic_f, #PETScContext::flag_hyperbolic_df,
  #PETScContext::flag_hyperbolic, #PETScContext::flag_parabolic, and #PETScContext::flag_source).

  This function computes \f${\bf F}\left({\bf U}\right)\f$, given \f${\bf U}\f$.

  \sa PetscIFunctionIMEX(), 
      TSSetRHSFunction (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSFunction.html)

  Note:
  + \f${\bf U}\f$ is \a Y in the code. PETsc denotes the state vector with \f${\bf Y}\f$ in its time integrators.
  + It is assumed that the reader is familiar with PETSc's implementation of IMEX time integrators, for
    example, TSARKIMEX (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX.html).
  + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html for documentation on
    PETSc's time integrators.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscRHSFunctionIMEX(
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
  IERR TransferVecFromPETSc(u,Y,context); CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t); CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                               solver->ghosts,mpi,u); CHECKERR(ierr);

  /* initialize right-hand side to zero */
  _ArraySetValue_(rhs,size*solver->nvars,0.0);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  if ((!strcmp(solver->SplitHyperbolicFlux,"yes")) && solver->flag_fdf_specified) {
    if (context->flag_hyperbolic_f == _EXPLICIT_) {
      IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->FdFFunction,solver->UpwindFdF);  
      CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    } 
    if (context->flag_hyperbolic_df == _EXPLICIT_) {
      IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF); 
      CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    }
  } else if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (context->flag_hyperbolic_f == _EXPLICIT_) {
      IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->FFunction,solver->Upwind);  
      CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
      IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF); 
      CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp, 1.0,rhs,size*solver->nvars);
    } 
    if (context->flag_hyperbolic_df == _EXPLICIT_) {
      IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF); 
      CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    }
  } else {
    if (context->flag_hyperbolic == _EXPLICIT_) {
      IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,solver->FFunction,solver->Upwind);  
      CHECKERR(ierr);
      _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    }
  }
  if (context->flag_parabolic == _EXPLICIT_) {
    IERR solver->ParabolicFunction (solver->par,u,solver,mpi,t);                        
    CHECKERR(ierr);
    _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
  }
  if (context->flag_source == _EXPLICIT_) {
    IERR solver->SourceFunction    (solver->source,u,solver,mpi,t);                     
    CHECKERR(ierr);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);
  }

  /* Transfer RHS to PETSc vector */
  IERR TransferVecToPETSc(rhs,F,context); CHECKERR(ierr);

  PetscFunctionReturn(0);
}

#endif
