/*! @file TimeRHSFunctionExplicit.c
    @brief Right-hand-side computation for explicit time integration
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*!
  This function computes the right-hand-side of the ODE given by
  \f{equation}{
    \frac {{\bf u}}{dt} = {\bf F}\left({\bf u}\right)
  \f}
  for explicit time integration methods, i.e., where
  \f{equation}{
    {\bf F}\left({\bf u}\right) = - {\bf F}_{\rm hyperbolic}\left({\bf u}\right) 
                                  + {\bf F}_{\rm parabolic} \left({\bf u}\right)
                                  + {\bf F}_{\rm source}    \left({\bf u}\right),
  \f}
  given the solution \f${\bf u}\f$ and the current simulation time.
*/
int TimeRHSFunctionExplicit(
                              double  *rhs, /*!< Array to hold the computed right-hand-side */
                              double  *u,   /*!< Array holding the solution */
                              void    *s,   /*!< Solver object of type #HyPar */
                              void    *m,   /*!< MPI object of type #MPIVariables */
                              double  t     /*!< Current simulation time */
                           ) 
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             d;
  _DECLARE_IERR_;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);                        CHECKERR(ierr);
  IERR solver->ApplyIBConditions(solver,mpi,u,t);                                   CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,u);                             CHECKERR(ierr);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->FFunction,solver->Upwind);
                                                                                    CHECKERR(ierr);
  IERR solver->ParabolicFunction (solver->par,u,solver,mpi,t);                      CHECKERR(ierr);
  IERR solver->SourceFunction    (solver->source,u,solver,mpi,t);                   CHECKERR(ierr);

  _ArraySetValue_(rhs,size*solver->nvars,0.0);
  _ArrayAXPY_(solver->hyp   ,-1.0,rhs,size*solver->nvars);
  _ArrayAXPY_(solver->par   , 1.0,rhs,size*solver->nvars);
  _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

  return(0);
}
