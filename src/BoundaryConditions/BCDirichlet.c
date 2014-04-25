#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCDirichletU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,size  ,indexb,boundary->is,ghosts,p);
      for (v=0; v<nvars; v++) phi[nvars*p+v] = boundary->DirichletValue[v];
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

int BCDirichletDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,size  ,indexb,boundary->is,ghosts,p);
      for (v=0; v<nvars; v++) phi[nvars*p+v] = 0.0;
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

