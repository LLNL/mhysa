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
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_], 
         DL[_MODEL_NVARS_*_MODEL_NVARS_], modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);

      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param); 

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_ (uavg,R,param,0);

      for (k=0; k<_MODEL_NVARS_; k++) D[k*_MODEL_NVARS_+k] = absolute(D[k*_MODEL_NVARS_+k]);
      MatMult(3,DL,D,L);
      MatMult(3,modA,R,DL);

      udiss[0] = modA[0*_MODEL_NVARS_+0]*udiff[0] + modA[0*_MODEL_NVARS_+1]*udiff[1] + modA[0*_MODEL_NVARS_+2]*udiff[2];
      udiss[1] = modA[1*_MODEL_NVARS_+0]*udiff[0] + modA[1*_MODEL_NVARS_+1]*udiff[1] + modA[1*_MODEL_NVARS_+2]*udiff[2];
      udiss[2] = modA[2*_MODEL_NVARS_+0]*udiff[0] + modA[2*_MODEL_NVARS_+1]*udiff[1] + modA[2*_MODEL_NVARS_+2]*udiff[2];

      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
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
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* calculate characteristic fluxes and variables */
      MatVecMult(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,0);
        eigL = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,0);
        eigR = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = D[k*_MODEL_NVARS_+k];

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
      MatVecMult(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
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
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Local Lax-Friedrich upwinding scheme */

      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* calculate characteristic fluxes and variables */
      MatVecMult(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,0);
        eigL = absolute(D[k*_MODEL_NVARS_+k]);
        _Euler1DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,0);
        eigR = absolute(D[k*_MODEL_NVARS_+k]);
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = absolute(D[k*_MODEL_NVARS_+k]);

        double alpha = max3(eigL,eigC,eigR);
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      IERR MatVecMult(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc); CHECKERR(ierr);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
