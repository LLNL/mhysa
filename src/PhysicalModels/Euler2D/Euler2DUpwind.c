#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

int Euler2DGetFlowVar        (double*,double*,double*,double*,double*,double*,void*);
int Euler2DRoeAverage        (double*,double*,double*,void*);
int Euler2DEigenvalues       (double*,double*,void*,int);
int Euler2DLeftEigenvectors  (double*,double*,void*,int);
int Euler2DRightEigenvectors (double*,double*,void*,int);

int Euler2DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      done,k,v;
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

      for (v=0; v<nvars; v++) udiff[v] = 0.5 * (uR[nvars*p+v] - uL[nvars*p+v]);

      IERR Euler2DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      IERR Euler2DEigenvalues       (uavg,D,param,dir); CHECKERR(ierr);
      IERR Euler2DLeftEigenvectors  (uavg,L,param,dir); CHECKERR(ierr);
      IERR Euler2DRightEigenvectors (uavg,R,param,dir); CHECKERR(ierr);

      for (k=0; k<nvars; k++) D[k*nvars+k] = absolute(D[k*nvars+k]);
      IERR MatMult(nvars,DL,D,L);    CHECKERR(ierr);
      IERR MatMult(nvars,modA,R,DL); CHECKERR(ierr);

      for (k=0; k<nvars; k++) {
        udiss[k] = 0; for (v=0; v<nvars; v++) udiss[k] += modA[k*nvars+v]*udiff[v];
      }
      
      for (v=0; v<nvars; v++) fI[nvars*p+v] = 0.5 * (fL[nvars*p+v]+fR[nvars*p+v]) - udiss[v];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler2DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      done,k;
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

      IERR Euler2DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      IERR Euler2DEigenvalues       (uavg,D,param,dir); CHECKERR(ierr);
      IERR Euler2DLeftEigenvectors  (uavg,L,param,dir); CHECKERR(ierr);
      IERR Euler2DRightEigenvectors (uavg,R,param,dir); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      IERR MatVecMult(nvars,ucL,L,&uL[nvars*p]);CHECKERR(ierr);
      IERR MatVecMult(nvars,ucR,L,&uR[nvars*p]);CHECKERR(ierr);
      IERR MatVecMult(nvars,fcL,L,&fL[nvars*p]);CHECKERR(ierr);
      IERR MatVecMult(nvars,fcR,L,&fR[nvars*p]);CHECKERR(ierr);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        IERR Euler2DEigenvalues(&uL[nvars*p],D,param,dir); CHECKERR(ierr);
        eigL = D[k*nvars+k];
        IERR Euler2DEigenvalues(&uR[nvars*p],D,param,dir); CHECKERR(ierr);
        eigR = D[k*nvars+k];
        IERR Euler2DEigenvalues(uavg        ,D,param,dir); CHECKERR(ierr);
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

int Euler2DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      done,k;
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

      /* Local Lax-Friedrich upwinding scheme */

      IERR Euler2DRoeAverage(uavg,&uL[nvars*p],&uR[nvars*p],param);  CHECKERR(ierr);

      IERR Euler2DEigenvalues       (uavg,D,param,dir); CHECKERR(ierr);
      IERR Euler2DLeftEigenvectors  (uavg,L,param,dir); CHECKERR(ierr);
      IERR Euler2DRightEigenvectors (uavg,R,param,dir); CHECKERR(ierr);

      /* calculate characteristic fluxes and variables */
      IERR MatVecMult(nvars,ucL,L,&uL[nvars*p]);CHECKERR(ierr);
      IERR MatVecMult(nvars,ucR,L,&uR[nvars*p]);CHECKERR(ierr);
      IERR MatVecMult(nvars,fcL,L,&fL[nvars*p]);CHECKERR(ierr);
      IERR MatVecMult(nvars,fcR,L,&fR[nvars*p]);CHECKERR(ierr);

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        IERR Euler2DEigenvalues(&uL[nvars*p],D,param,dir); CHECKERR(ierr);
        eigL = D[k*nvars+k];
        IERR Euler2DEigenvalues(&uR[nvars*p],D,param,dir); CHECKERR(ierr);
        eigR = D[k*nvars+k];
        IERR Euler2DEigenvalues(uavg        ,D,param,dir); CHECKERR(ierr);
        eigC = D[k*nvars+k];

        double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      IERR MatVecMult(nvars,&fI[nvars*p],R,fc); CHECKERR(ierr);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
