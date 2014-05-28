#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int TimeRHSFunctionExplicit(double *rhs,double *u,void *s,void *m, double t) 
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             d;
  _DECLARE_IERR_;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);          CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                 solver->ghosts,mpi,u);                 CHECKERR(ierr);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
  _ArraySetValue_(rhs,size*solver->nvars,0.0);
  if (solver->HyperbolicFunction) {
    IERR solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,NULL);   CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp    ,-1.0,rhs,size*solver->nvars);
  }
  if (solver->HyperbolicFunction1) {
    IERR solver->HyperbolicFunction1(solver->hyp,u,solver,mpi,t,NULL);  CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp    ,-1.0,rhs,size*solver->nvars);
  }
  if (solver->HyperbolicFunction2) {
    IERR solver->HyperbolicFunction2(solver->hyp,u,solver,mpi,t,NULL);  CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp    ,-1.0,rhs,size*solver->nvars);
  }
  if (solver->ParabolicFunction) {
    IERR solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    _ArrayAXPY_(solver->par    , 1.0,rhs,size*solver->nvars);
  }
  if (solver->SourceFunction) {
    IERR solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    _ArrayAXPY_(solver->source , 1.0,rhs,size*solver->nvars);
  }

  return(0);
}
