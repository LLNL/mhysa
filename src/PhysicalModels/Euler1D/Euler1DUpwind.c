#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DGetFlowVar        (double*,double*,double*,double*,double*,void*);
int Euler1DRoeAverage        (double*,double*,double*,void*);
int Euler1DEigenvalues       (double*,double*,void*,int);
int Euler1DLeftEigenvectors  (double*,double*,void*,int);
int Euler1DRightEigenvectors (double*,double*,void*,int);

int Euler1DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars], DL[nvars*nvars], modA[nvars*nvars];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double udiff[nvars], uavg[nvars],udiss[nvars];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[nvars*p+0] - uL[nvars*p+0]);
      udiff[1] = 0.5 * (uR[nvars*p+1] - uL[nvars*p+1]);
      udiff[2] = 0.5 * (uR[nvars*p+2] - uL[nvars*p+2]);

      IERR Euler1DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param); 

      IERR Euler1DEigenvalues       (uavg,D,param,0); CHECKERR(ierr);
      IERR Euler1DLeftEigenvectors  (uavg,L,param,0); CHECKERR(ierr);
      IERR Euler1DRightEigenvectors (uavg,R,param,0); CHECKERR(ierr);

      for (k=0; k<nvars; k++) D[k*nvars+k] = absolute(D[k*nvars+k]);
      IERR MatMult(3,DL,D,L);    CHECKERR(ierr);
      IERR MatMult(3,modA,R,DL); CHECKERR(ierr);

      udiss[0] = modA[0*nvars+0]*udiff[0] + modA[0*nvars+1]*udiff[1] + modA[0*nvars+2]*udiff[2];
      udiss[1] = modA[1*nvars+0]*udiff[0] + modA[1*nvars+1]*udiff[1] + modA[1*nvars+2]*udiff[2];
      udiss[2] = modA[2*nvars+0]*udiff[0] + modA[2*nvars+1]*udiff[1] + modA[2*nvars+2]*udiff[2];

      fI[nvars*p+0] = 0.5 * (fL[nvars*p+0]+fR[nvars*p+0]) - udiss[0];
      fI[nvars*p+1] = 0.5 * (fL[nvars*p+1]+fR[nvars*p+1]) - udiss[1];
      fI[nvars*p+2] = 0.5 * (fL[nvars*p+2]+fR[nvars*p+2]) - udiss[2];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler1DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double uavg[nvars], fcL[nvars], fcR[nvars], ucL[nvars], ucR[nvars], fc[nvars];

      /* Roe-Fixed upwinding scheme */

      IERR Euler1DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      IERR Euler1DEigenvalues       (uavg,D,param,0); CHECKERR(ierr);
      IERR Euler1DLeftEigenvectors  (uavg,L,param,0); CHECKERR(ierr);
      IERR Euler1DRightEigenvectors (uavg,R,param,0); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      IERR MatVecMult(nvars,ucL,L,&uL[nvars*p]); CHECKERR(ierr);
      IERR MatVecMult(nvars,ucR,L,&uR[nvars*p]); CHECKERR(ierr);
      IERR MatVecMult(nvars,fcL,L,&fL[nvars*p]); CHECKERR(ierr);
      IERR MatVecMult(nvars,fcR,L,&fR[nvars*p]); CHECKERR(ierr);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        IERR Euler1DEigenvalues(&uL[nvars*p],D,param,0); CHECKERR(ierr);
        eigL = D[k*nvars+k];
        IERR Euler1DEigenvalues(&uR[nvars*p],D,param,0); CHECKERR(ierr);
        eigR = D[k*nvars+k];
        IERR Euler1DEigenvalues(uavg        ,D,param,0); CHECKERR(ierr);
        eigC = D[k*nvars+k];

        if ((eigL > 0) && (eigC > 0) && (eigR > 0)) {
          fc[k] = fcL[k];
        } else if ((eigL < 0) && (eigC < 0) && (eigR < 0)) {
          fc[k] = fcR[k];
        } else {
          double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      IERR MatVecMult(nvars,&fI[nvars*p],R,fc); CHECKERR(ierr);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler1DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double uavg[nvars], fcL[nvars], fcR[nvars], ucL[nvars], ucR[nvars], fc[nvars];

      /* Roe-Fixed upwinding scheme */

      IERR Euler1DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      IERR Euler1DEigenvalues       (uavg,D,param,0); CHECKERR(ierr);
      IERR Euler1DLeftEigenvectors  (uavg,L,param,0); CHECKERR(ierr);
      IERR Euler1DRightEigenvectors (uavg,R,param,0); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      IERR MatVecMult(nvars,ucL,L,&uL[nvars*p]); CHECKERR(ierr);
      IERR MatVecMult(nvars,ucR,L,&uR[nvars*p]); CHECKERR(ierr);
      IERR MatVecMult(nvars,fcL,L,&fL[nvars*p]); CHECKERR(ierr);
      IERR MatVecMult(nvars,fcR,L,&fR[nvars*p]); CHECKERR(ierr);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        IERR Euler1DEigenvalues(&uL[nvars*p],D,param,0); CHECKERR(ierr);
        eigL = absolute(D[k*nvars+k]);
        IERR Euler1DEigenvalues(&uR[nvars*p],D,param,0); CHECKERR(ierr);
        eigR = absolute(D[k*nvars+k]);
        IERR Euler1DEigenvalues(uavg        ,D,param,0); CHECKERR(ierr);
        eigC = absolute(D[k*nvars+k]);

        double alpha = max3(eigL,eigC,eigR);
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      IERR MatVecMult(nvars,&fI[nvars*p],R,fc); CHECKERR(ierr);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
