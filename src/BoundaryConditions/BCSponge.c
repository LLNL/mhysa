#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

int BCSpongeSource(void *b,int ndims,int nvars,int ghosts,int *size,double *grid,double *u,double *source)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int            dim       = boundary->dim;
  int            face      = boundary->face;
  double         *uref     = boundary->SpongeValue;
  double         *xmin     = boundary->xmin;
  double         *xmax     = boundary->xmax;
  int            v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int i, istart, iend;
      i       = indexb[dim];
      istart  = boundary->is[dim];
      iend    = boundary->ie[dim];
      double x, xstart, xend;
      _GetCoordinate_(dim,i,size,ghosts,grid,x);
      xstart = xmin[dim];
      xend   = xmax[dim];
      /* calculate sigma */
      double sigma;
      if (face > 0) sigma = (x - xstart) / (xend - xstart);
      else          sigma = (x - xend  ) / (xstart - xend);
      /* add to the source term */
      int p; _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p);
      for (v=0; v<nvars; v++) 
        source[nvars*p+v] -= (sigma * (u[nvars*p+v]-uref[nvars*p+v]));
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

/*
 * These are dummy functions that don't do anything. They are called by ApplyBoundaryConditions()
 * when going through all the boundary zones. The actual implementation of the sponge BC is through
 * the source function that calls the function above.
*/ 

int BCSpongeUDummy(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double waqt)
{
  return(0);
}

int BCSpongeDUDummy(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref,double waqt)
{
  return(0);
}
