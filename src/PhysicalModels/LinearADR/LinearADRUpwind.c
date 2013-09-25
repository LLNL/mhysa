#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <hypar.h>

int LinearADRUpwind(double *fI,double *fL,double *fR,double *u,int dir,void *s)
{
  HyPar     *solver = (HyPar*)      s;
  LinearADR *param  = (LinearADR*)  solver->physics;
  int       ierr    = 0,done,v;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  double *a = param->a;

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
      for (v = 0; v < nvars; v++)  
        fI[nvars*p+v] = (a[nvars*dir+v] > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  free(index_inter);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);

  return(0);
}
