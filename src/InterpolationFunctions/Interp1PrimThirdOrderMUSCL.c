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

int Interp1PrimThirdOrderMUSCL(double *fI,double *fC,double *u,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  MUSCLParameters *muscl  = (MUSCLParameters*) solver->interp;
  int             ierr    = 0;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  double one_third = 1.0/3.0;
  double one_sixth = 1.0/6.0;

  if ((!fI) || (!fC)) {
    fprintf(stderr, "Error in Interp1PrimThirdOrderMUSCL(): input arrays not allocated.\n");
    return(1);
  }
  if (ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in Interp1PrimThirdOrderMUSCL(): insufficient number of ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  int done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  if (upw > 0) {
    while (!done) {
      ierr = ArrayCopy1D_int(index_outer,indexC,ndims); CHECKERR(ierr);
      ierr = ArrayCopy1D_int(index_outer,indexI,ndims); CHECKERR(ierr);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p = ArrayIndex1D(ndims,bounds_inter,indexI,NULL,0);
        int qm1,qm2,qp1,v;
        indexC[dir] = indexI[dir]-2; qm2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-1; qm1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]  ; qp1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
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
      done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
    }
  } else {
    while (!done) {
      ierr = ArrayCopy1D_int(index_outer,indexC,ndims); CHECKERR(ierr);
      ierr = ArrayCopy1D_int(index_outer,indexI,ndims); CHECKERR(ierr);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p = ArrayIndex1D(ndims,bounds_inter,indexI,NULL,0);
        int qm1,qp1,qp2,v;
        indexC[dir] = indexI[dir]-1; qm1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]  ; qp1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]+1; qp2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
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
      done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
    }
  }
  
  return(0);
}
