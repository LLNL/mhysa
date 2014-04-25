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

int BCPeriodicU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

  if ((boundary->on_this_proc) && (mpi->iproc[dim] == 1)) {
    int bounds[ndims], index1[ndims], index2[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(index1,ndims,0);
    _ArraySetValue_(index2,ndims,0);
    int done = 0;
    while (!done) {
      int p1 = 0, p2 = 0;
      _ArrayCopy1D_(index1,index2,ndims);
      if (face == 1) {
        index2[dim] = index1[dim] + size[dim]-ghosts;
        _ArrayIndex1DWO_(ndims,size,index1,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,index2,ghosts,p2);
      } else if (face == -1) {
        _ArrayIndex1DWO_(ndims,size,index1,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,index1,ghosts,p2);
      }
      for (v=0; v<nvars; v++) phi[nvars*p1+v] = phi[nvars*p2+v];
      _ArrayIncrementIndex_(ndims,bounds,index1,done);
    }
  }
  return(0);
}

int BCPeriodicDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

  if ((boundary->on_this_proc) && (mpi->iproc[dim] == 1)) {
    int bounds[ndims], index1[ndims], index2[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(index1,ndims,0);
    _ArraySetValue_(index2,ndims,0);
    int done = 0;
    while (!done) {
      int p1 = 0, p2 = 0;
      _ArrayCopy1D_(index1,index2,ndims);
      if (face == 1) {
        index2[dim] = index1[dim] + size[dim]-ghosts;
        _ArrayIndex1DWO_(ndims,size,index1,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,index2,ghosts,p2);
      } else if (face == -1) {
        _ArrayIndex1DWO_(ndims,size,index1,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,index1,ghosts,p2);
      }
      for (v=0; v<nvars; v++) phi[nvars*p1+v] = phi[nvars*p2+v];
      _ArrayIncrementIndex_(ndims,bounds,index1,done);
    }
  }
  return(0);
}
