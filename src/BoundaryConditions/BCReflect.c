#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCReflectU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int var   = boundary->var;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0);
    int done = 0;
    while (!done) {
      int p1, p2;
      _ArrayCopy1D_(indexb,indexi,ndims);
      _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
      if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
      else return(1);
      _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
      phi[nvars*p1+var] = -phi[nvars*p2+var];
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

int BCReflectDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int var   = boundary->var;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0);
    int done = 0;
    while (!done) {
      int p1, p2;
      _ArrayCopy1D_(indexb,indexi,ndims);
      _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
      if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
      else return(1);
      _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
      phi[nvars*p1+var] = -phi[nvars*p2+var];
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}
