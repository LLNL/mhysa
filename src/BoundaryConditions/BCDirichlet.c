#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCDirichlet(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            var       = boundary->var;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,size  ,indexb,boundary->is,ghosts,p);
      phi[nvars*p+var] = boundary->DirichletValue[var];
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

