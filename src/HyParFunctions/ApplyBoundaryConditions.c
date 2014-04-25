#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>

/*
 * The argument flag:
 * 1 -> argument x is (\Delta U), hence apply the BCs for
 *      (\Delta U)
 * 0 -> argument x is U, hence apply the BCs for U
*/

int ApplyBoundaryConditions(void *s,void *m,double *x,double *xref,int flag,double waqt) 
{
  HyPar           *solver   = (HyPar*)          s;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  int             nb        = solver->nBoundaryZones;
  _DECLARE_IERR_;


  /* Apply domain boundary conditions to p */
  int n;
  if (flag) {
    for (n = 0; n < nb; n++) {
      IERR boundary[n].BCFunctionDU(&boundary[n],mpi,solver->ndims,solver->nvars,
                                   solver->dim_local,solver->ghosts,x,xref,waqt);
      CHECKERR(ierr);
    }
  } else {
    for (n = 0; n < nb; n++) {
      IERR boundary[n].BCFunctionU(&boundary[n],mpi,solver->ndims,solver->nvars,
                                   solver->dim_local,solver->ghosts,x,waqt);
      CHECKERR(ierr);
    }
  }

  return(0);
}
