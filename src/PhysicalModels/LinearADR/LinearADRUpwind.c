#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <hypar.h>

int LinearADRUpwind(double *fI,double *fL,double *fR,double *uL,double *uR,
                    double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)      s;
  LinearADR *param  = (LinearADR*)  solver->physics;
  int       done,v;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  double *a = param->a;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (v = 0; v < nvars; v++)  
        fI[nvars*p+v] = (a[nvars*dir+v] > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
