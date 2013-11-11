#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

inline int Euler1DGetFlowVar        (double*,double*,double*,double*,double*,void*);
inline int Euler1DRoeAverage        (double*,double*,double*,void*);
inline int Euler1DEigenvalues       (double*,double*,void*,int);
inline int Euler1DLeftEigenvectors  (double*,double*,void*,int);
inline int Euler1DRightEigenvectors (double*,double*,void*,int);

int Euler1DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       ierr    = 0,done,k;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars], DL[nvars*nvars], modA[nvars*nvars];

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double udiff[nvars], uavg[nvars],udiss[nvars];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[nvars*p+0] - uL[nvars*p+0]);
      udiff[1] = 0.5 * (uR[nvars*p+1] - uL[nvars*p+1]);
      udiff[2] = 0.5 * (uR[nvars*p+2] - uL[nvars*p+2]);

      ierr = Euler1DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      ierr = Euler1DEigenvalues       (uavg,D,param,0); CHECKERR(ierr);
      ierr = Euler1DLeftEigenvectors  (uavg,L,param,0); CHECKERR(ierr);
      ierr = Euler1DRightEigenvectors (uavg,R,param,0); CHECKERR(ierr);

      for (k=0; k<nvars; k++) D[k*nvars+k] = absolute(D[k*nvars+k]);
      ierr = MatMult(3,DL,D,L);    CHECKERR(ierr);
      ierr = MatMult(3,modA,R,DL); CHECKERR(ierr);

      udiss[0] = modA[0*nvars+0]*udiff[0] + modA[0*nvars+1]*udiff[1] + modA[0*nvars+2]*udiff[2];
      udiss[1] = modA[1*nvars+0]*udiff[0] + modA[1*nvars+1]*udiff[1] + modA[1*nvars+2]*udiff[2];
      udiss[2] = modA[2*nvars+0]*udiff[0] + modA[2*nvars+1]*udiff[1] + modA[2*nvars+2]*udiff[2];

      fI[nvars*p+0] = 0.5 * (fL[nvars*p+0]+fR[nvars*p+0]) - udiss[0];
      fI[nvars*p+1] = 0.5 * (fL[nvars*p+1]+fR[nvars*p+1]) - udiss[1];
      fI[nvars*p+2] = 0.5 * (fL[nvars*p+2]+fR[nvars*p+2]) - udiss[2];
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  return(0);
}

int Euler1DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       ierr    = 0,done,k;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars];

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double uavg[nvars], fcL[nvars], fcR[nvars], ucL[nvars], ucR[nvars], fc[nvars];

      /* Roe-Fixed upwinding scheme */

      ierr = Euler1DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      ierr = Euler1DEigenvalues       (uavg,D,param,0); CHECKERR(ierr);
      ierr = Euler1DLeftEigenvectors  (uavg,L,param,0); CHECKERR(ierr);
      ierr = Euler1DRightEigenvectors (uavg,R,param,0); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      ierr = MatVecMult(nvars,ucL,L,&uL[nvars*p]);
      ierr = MatVecMult(nvars,ucR,L,&uR[nvars*p]);
      ierr = MatVecMult(nvars,fcL,L,&fL[nvars*p]);
      ierr = MatVecMult(nvars,fcR,L,&fR[nvars*p]);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        ierr = Euler1DEigenvalues(&uL[nvars*p],D,param,0); CHECKERR(ierr);
        eigL = D[k*nvars+k];
        ierr = Euler1DEigenvalues(&uR[nvars*p],D,param,0); CHECKERR(ierr);
        eigR = D[k*nvars+k];
        ierr = Euler1DEigenvalues(uavg        ,D,param,0); CHECKERR(ierr);
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
      ierr = MatVecMult(nvars,&fI[nvars*p],R,fc); CHECKERR(ierr);
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  return(0);
}

int Euler1DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       ierr    = 0,done,k;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars];

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double uavg[nvars], fcL[nvars], fcR[nvars], ucL[nvars], ucR[nvars], fc[nvars];

      /* Roe-Fixed upwinding scheme */

      ierr = Euler1DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      ierr = Euler1DEigenvalues       (uavg,D,param,0); CHECKERR(ierr);
      ierr = Euler1DLeftEigenvectors  (uavg,L,param,0); CHECKERR(ierr);
      ierr = Euler1DRightEigenvectors (uavg,R,param,0); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      ierr = MatVecMult(nvars,ucL,L,&uL[nvars*p]);
      ierr = MatVecMult(nvars,ucR,L,&uR[nvars*p]);
      ierr = MatVecMult(nvars,fcL,L,&fL[nvars*p]);
      ierr = MatVecMult(nvars,fcR,L,&fR[nvars*p]);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        ierr = Euler1DEigenvalues(&uL[nvars*p],D,param,0); CHECKERR(ierr);
        eigL = absolute(D[k*nvars+k]);
        ierr = Euler1DEigenvalues(&uR[nvars*p],D,param,0); CHECKERR(ierr);
        eigR = absolute(D[k*nvars+k]);
        ierr = Euler1DEigenvalues(uavg        ,D,param,0); CHECKERR(ierr);
        eigC = absolute(D[k*nvars+k]);

        double alpha = max3(eigL,eigC,eigR);
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      ierr = MatVecMult(nvars,&fI[nvars*p],R,fc); CHECKERR(ierr);
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  return(0);
}
