#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_omp
#include <omp.h>
#endif

/* 
  Fifth order WENO interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int Interp1PrimFifthOrderWENO(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  double *ww1, *ww2, *ww3;
  ww1 = weno->w1 + (upw < 0 ? 2*weno->size : 0) + (u==fC ? weno->size : 0) + weno->offset[dir];
  ww2 = weno->w2 + (upw < 0 ? 2*weno->size : 0) + (u==fC ? weno->size : 0) + weno->offset[dir];
  ww3 = weno->w3 + (upw < 0 ? 2*weno->size : 0) + (u==fC ? weno->size : 0) + weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  int i;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }
      int v; 
      for (v=0; v<nvars; v++)  {
        /* Defining stencil points */
        double fm3, fm2, fm1, fp1, fp2;
        fm3 = fC[qm3*nvars+v];
        fm2 = fC[qm2*nvars+v];
        fm1 = fC[qm1*nvars+v];
        fp1 = fC[qp1*nvars+v];
        fp2 = fC[qp2*nvars+v];

        /* Candidate stencils and their optimal weights*/
        double f1, f2, f3;
        f1 = (2*one_sixth)*fm3 - (7.0*one_sixth)*fm2 + (11.0*one_sixth)*fm1;
        f2 = (-one_sixth)*fm2 + (5.0*one_sixth)*fm1 + (2*one_sixth)*fp1;
        f3 = (2*one_sixth)*fm1 + (5*one_sixth)*fp1 - (one_sixth)*fp2;

        /* calculate WENO weights */
        double w1,w2,w3;
        w1 = *(ww1+p*nvars+v);
        w2 = *(ww2+p*nvars+v);
        w3 = *(ww3+p*nvars+v);

        fI[p*nvars+v] = w1*f1 + w2*f2 + w3*f3;
      }
    }
  }

  return(0);
}
