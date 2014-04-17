#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>

int ApplyBoundaryConditions(void *s,void *m,double *x) 
{
  HyPar           *solver   = (HyPar*)          s;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  int             nb        = solver->nBoundaryZones;
  _DECLARE_IERR_;


  /* Apply domain boundary conditions to p */
  int n;
  for (n = 0; n < nb; n++) {
    IERR boundary[n].BCFunctionU(&boundary[n],mpi,solver->ndims,solver->nvars,
                                  solver->dim_local,solver->ghosts,x);
    CHECKERR(ierr);
  }

  return(0);
}
