#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fpdoublewell.h>
#include <hypar.h>

int FPDoubleWellUpwind(double *fI,double *fL,double *fR,double *u,int dir,void *s)
{
  HyPar        *solver = (HyPar*) s;
  int          ierr    = 0,done,v;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ghosts  = solver->ghosts;
  int *dim    = solver->dim_local;

  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *index_inter  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  int *bounds_inter = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double x = 0.5 * (solver->x[index_inter[dir]-1+ghosts] + solver->x[index_inter[dir]+ghosts]);
      for (v = 0; v < nvars; v++)  
        fI[nvars*p+v] = (drift(x) > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  free(index_inter);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);

  return(0);
}
