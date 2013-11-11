#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCDirichlet(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            ierr      = 0;

  int var   = boundary->var;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    ierr = ArraySubtract1D_int(bounds,boundary->ie,boundary->is,ndims); CHECKERR(ierr);
    ierr = ArraySetValue_int  (indexb,ndims,0);                         CHECKERR(ierr);
    int done = 0;
    while (!done) {
      int p = ArrayIndex1D(ndims,size  ,indexb,boundary->is,ghosts);
      phi[nvars*p+var] = boundary->DirichletValue[var];
      done = ArrayIncrementIndex(ndims,bounds,indexb);
    }
  }
  return(0);
}

