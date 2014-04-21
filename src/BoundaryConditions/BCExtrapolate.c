#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCExtrapolateU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0);
    int done = 0;
    while (!done) {
      _ArrayCopy1D_(indexb,indexi,ndims);
      _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
      if (face == 1)        indexi[dim] = 0;
      else if (face == -1)  indexi[dim] = size[dim]-1;
      else                  return(1);
      int p1,p2;
      _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
      for (v=0; v<nvars; v++) phi[nvars*p1+v] = phi[nvars*p2+v];
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

int BCExtrapolateDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0);
    int done = 0;
    while (!done) {
      _ArrayCopy1D_(indexb,indexi,ndims);
      _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
      if (face == 1)        indexi[dim] = 0;
      else if (face == -1)  indexi[dim] = size[dim]-1;
      else                  return(1);
      int p1,p2;
      _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
      for (v=0; v<nvars; v++) phi[nvars*p1+v] = phi[nvars*p2+v];
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}
