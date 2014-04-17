#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Third order MUSCL interpolation with Koren's limiter (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 2

int Interp1PrimThirdOrderMUSCL(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  MUSCLParameters *muscl  = (MUSCLParameters*) solver->interp;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  double one_third = 1.0/3.0;
  double one_sixth = 1.0/6.0;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  if (upw > 0) {
    while (!done) {
      _ArrayCopy1D_(index_outer,indexC,ndims);
      _ArrayCopy1D_(index_outer,indexI,ndims);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
        int qm1,qm2,qp1,v;
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        for (v=0; v<nvars; v++)  {
          /* Defining stencil points */
          double m2, m1, p1;
          m2 = fC[qm2*nvars+v];
          m1 = fC[qm1*nvars+v];
          p1 = fC[qp1*nvars+v];

          double fdiff = p1 - m1;
          double bdiff = m1 - m2;
          double limit =  (3*fdiff*bdiff + muscl->eps) 
                        / (2*(fdiff-bdiff)*(fdiff-bdiff) + 3*fdiff*bdiff + muscl->eps);

          fI[p*nvars+v] = m1 +  limit * (one_third*fdiff + one_sixth*bdiff);
        }
      }
      _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
    }
  } else {
    while (!done) {
      _ArrayCopy1D_(index_outer,indexC,ndims);
      _ArrayCopy1D_(index_outer,indexI,ndims);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
        int qm1,qp1,qp2,v;
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
        for (v=0; v<nvars; v++)  {
          /* Defining stencil points */
          double m1, p1, p2;
          m1 = fC[qm1*nvars+v];
          p1 = fC[qp1*nvars+v];
          p2 = fC[qp2*nvars+v];

          double fdiff = p2 - p1;
          double bdiff = p1 - m1;
          double limit =  (3*fdiff*bdiff + muscl->eps) 
                        / (2*(fdiff-bdiff)*(fdiff-bdiff) + 3*fdiff*bdiff + muscl->eps);

          fI[p*nvars+v] = p1 -  limit * (one_third*fdiff + one_sixth*bdiff);
        }
      }
      _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
    }
  }
  
  return(0);
}
