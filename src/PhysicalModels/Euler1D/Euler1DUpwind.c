#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

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

      _Euler1DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param); 

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_ (uavg,R,param,0);

      for (k=0; k<nvars; k++) D[k*nvars+k] = absolute(D[k*nvars+k]);
      MatMult(3,DL,D,L);
      MatMult(3,modA,R,DL);

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

      _Euler1DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* calculate characteristic fluxes and variables */
      MatVecMult(nvars,ucL,L,(uL+nvars*p));
      MatVecMult(nvars,ucR,L,(uR+nvars*p));
      MatVecMult(nvars,fcL,L,(fL+nvars*p));
      MatVecMult(nvars,fcR,L,(fR+nvars*p));

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((uL+nvars*p),D,param,0);
        eigL = D[k*nvars+k];
        _Euler1DEigenvalues_((uR+nvars*p),D,param,0);
        eigR = D[k*nvars+k];
        _Euler1DEigenvalues_(uavg,D,param,0);
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
      MatVecMult(nvars,(fI+nvars*p),R,fc);
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

      /* Local Lax-Friedrich upwinding scheme */

      _Euler1DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* calculate characteristic fluxes and variables */
      MatVecMult(nvars,ucL,L,(uL+nvars*p));
      MatVecMult(nvars,ucR,L,(uR+nvars*p));
      MatVecMult(nvars,fcL,L,(fL+nvars*p));
      MatVecMult(nvars,fcR,L,(fR+nvars*p));

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((uL+nvars*p),D,param,0);
        eigL = absolute(D[k*nvars+k]);
        _Euler1DEigenvalues_((uR+nvars*p),D,param,0);
        eigR = absolute(D[k*nvars+k]);
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = absolute(D[k*nvars+k]);

        double alpha = max3(eigL,eigC,eigR);
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      IERR MatVecMult(nvars,(fI+nvars*p),R,fc); CHECKERR(ierr);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
