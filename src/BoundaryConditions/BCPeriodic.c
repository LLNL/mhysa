#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>

int BCPeriodic(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;
  int            ierr     = 0, d;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int var   = boundary->var;

  /* allocate buffer to fetch data from other side of the domain */
  int buffer_size = 1;
  for (d = 0; d < ndims; d++) {
    if (d == dim) buffer_size *= ghosts;
    else          buffer_size *= size[d];
  }
  double *buf = (double*) calloc (buffer_size*nvars,sizeof(double));

  /* first and last processes along the specified dimension participate */
  if ((mpi->ip[dim] == 0) || (mpi->ip[dim] == mpi->iproc[dim]-1)) {
    int *source, *dest, *limits;
    source = (int*) calloc (ndims  ,sizeof(int));
    dest   = (int*) calloc (ndims  ,sizeof(int));
    limits = (int*) calloc (2*ndims,sizeof(int));
    ierr = ArrayCopy1D_int(mpi->ip,source,ndims); CHECKERR(ierr);
    ierr = ArrayCopy1D_int(mpi->ip,dest  ,ndims); CHECKERR(ierr);

    if (face == 1) {
      source[dim] = mpi->iproc[dim]-1;
      dest  [dim] = 0;
      for (d=0; d<ndims; d++) {
        if (d == dim) {
          limits[2*d]   = size[dim]-ghosts;
          limits[2*d+1] = size[dim];
        } else {
          limits[2*d]   = 0;
          limits[2*d+1] = size[dim];
        }
      }
    } else if (face == -1) {
      source[dim] = 0;
      dest  [dim] = mpi->iproc[dim]-1;
      for (d=0; d<ndims; d++) {
        if (d == dim) {
          limits[2*d]   = 0;
          limits[2*d+1] = ghosts;
        } else {
          limits[2*d]   = 0;
          limits[2*d+1] = size[dim];
        }
      }
    }
    ierr = MPIGetArrayDatanD(buf,phi,source,dest,limits,size,ghosts,ndims,nvars,mpi); CHECKERR(ierr);
    free(source);
    free(dest);
    free(limits);
  }
  if (boundary->on_this_proc) {
    int *bounds = (int*) calloc (ndims,sizeof(int));
    int *index  = (int*) calloc (ndims,sizeof(int));
    ierr = ArraySubtract1D_int(bounds,boundary->ie,boundary->is,ndims); CHECKERR(ierr);
    ierr = ArraySetValue_int  (index,ndims,0);                          CHECKERR(ierr);
    int done = 0;
    while (!done) {
      int p1 = ArrayIndex1D(ndims,size  ,index,boundary->is,ghosts);
      int p2 = ArrayIndex1D(ndims,bounds,index,NULL        ,0     );
      phi[nvars*p1+var] = buf[nvars*p2+var];
      done = ArrayIncrementIndex(ndims,bounds,index);
    }
  }
  free(buf);
  return(0);
}

