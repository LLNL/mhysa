#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>

/*
  This function implements periodic BCs through ghost points
  only if the number of processes along that particular dimension
  is 1.

  If it's greater than 1, it's handled by MPIExchangeBoundaries()
*/

int BCPeriodic(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;
  int            ierr      = 0;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int var   = boundary->var;

  if ((boundary->on_this_proc) && (mpi->iproc[dim] == 1)) {
    int bounds[ndims], index1[ndims], index2[ndims];
    ierr = ArraySubtract1D_int(bounds,boundary->ie,boundary->is,ndims); CHECKERR(ierr);
    ierr = ArraySetValue_int  (index1,ndims,0);                         CHECKERR(ierr);
    ierr = ArraySetValue_int  (index2,ndims,0);                         CHECKERR(ierr);
    int done = 0;
    while (!done) {
      int p1 = 0, p2 = 0;
      ierr = ArrayCopy1D_int(index1,index2,ndims);CHECKERR(ierr);
      if (face == 1) {
        index2[dim] = index1[dim] + size[dim]-ghosts;
        p1 = ArrayIndex1D(ndims,size,index1,boundary->is,ghosts);
        p2 = ArrayIndex1D(ndims,size,index2,NULL        ,ghosts);
      } else if (face == -1) {
        p1 = ArrayIndex1D(ndims,size,index1,boundary->is,ghosts);
        p2 = ArrayIndex1D(ndims,size,index1,NULL        ,ghosts);
      }
      phi[nvars*p1+var] = phi[nvars*p2+var];
      done = ArrayIncrementIndex(ndims,bounds,index1);
    }
  }
  return(0);
}

