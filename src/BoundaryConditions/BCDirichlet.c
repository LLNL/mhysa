#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCDirichlet(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            ierr      = 0, d;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int var   = boundary->var;

  if (boundary->on_this_proc) {
    int *bounds = (int*) calloc (ndims,sizeof(int));
    int *indexb = (int*) calloc (ndims,sizeof(int));  /* boundary index */
    int *indexi = (int*) calloc (ndims,sizeof(int));  /* interior index */
    ierr = ArraySubtract1D_int(bounds,boundary->ie,boundary->is,ndims); CHECKERR(ierr);
    ierr = ArraySetValue_int  (indexb,ndims,0);                         CHECKERR(ierr);
    int done = 0;
    while (!done) {
      int p = ArrayIndex1D(ndims,size  ,indexb,boundary->is,ghosts);
      phi[nvars*p+var] = boundary->DirichletValue[var];
      done = ArrayIncrementIndex(ndims,bounds,indexb);
    }
    free(bounds);
    free(indexb);
    free(indexi);
  }
  return(0);
}

