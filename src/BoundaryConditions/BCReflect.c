#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCReflect(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            ierr      = 0;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int var   = boundary->var;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims], indexi[ndims];
    ierr = ArraySubtract1D_int(bounds,boundary->ie,boundary->is,ndims); CHECKERR(ierr);
    ierr = ArraySetValue_int  (indexb,ndims,0);                         CHECKERR(ierr);
    int done = 0;
    while (!done) {
      int p1, p2;
      ierr = ArrayCopy1D_int(indexb,indexi,ndims);              CHECKERR(ierr);
      ierr = ArrayAdd1D_int (indexi,indexi,boundary->is,ndims); CHECKERR(ierr);
      if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = size[dim]-indexb[dim];
      else return(1);
      p1 = ArrayIndex1D(ndims,size  ,indexb,boundary->is,ghosts);
      p2 = ArrayIndex1D(ndims,size  ,indexi,NULL        ,ghosts);
      phi[nvars*p1+var] = -phi[nvars*p2+var];
      done = ArrayIncrementIndex(ndims,bounds,indexb);
    }
  }
  return(0);
}

