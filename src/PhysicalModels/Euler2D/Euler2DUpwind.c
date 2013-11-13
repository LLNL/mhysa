#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

int Euler2DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      done,k,v;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
         L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
         modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      for (v=0; v<_MODEL_NVARS_; v++) udiff[v] = 0.5 * (uR[_MODEL_NVARS_*p+v] - uL[_MODEL_NVARS_*p+v]);

      _Euler2DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler2DEigenvalues_(uavg,D,param,dir);
      _Euler2DLeftEigenvectors_(uavg,L,param,dir);
      _Euler2DRightEigenvectors_(uavg,R,param,dir);

      for (k=0; k<_MODEL_NVARS_; k++) D[k*_MODEL_NVARS_+k] = absolute(D[k*_MODEL_NVARS_+k]);
      MatMult(_MODEL_NVARS_,DL,D,L);
      MatMult(_MODEL_NVARS_,modA,R,DL);

      for (k=0; k<_MODEL_NVARS_; k++) {
        udiss[k] = 0; for (v=0; v<_MODEL_NVARS_; v++) udiss[k] += modA[k*_MODEL_NVARS_+v]*udiff[v];
      }
      
      for (v=0; v<_MODEL_NVARS_; v++) 
        fI[_MODEL_NVARS_*p+v] = 0.5 * (fL[_MODEL_NVARS_*p+v]+fR[_MODEL_NVARS_*p+v]) - udiss[v];
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

      _Euler2DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler2DEigenvalues_(uavg,D,param,dir);
      _Euler2DLeftEigenvectors_ (uavg,L,param,dir);
      _Euler2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler2DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,dir);
        eigL = D[k*_MODEL_NVARS_+k];
        _Euler2DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,dir);
        eigR = D[k*_MODEL_NVARS_+k];
        _Euler2DEigenvalues_(uavg        ,D,param,dir);
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

int Euler2DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      done,k;

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

      _Euler2DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler2DEigenvalues_(uavg,D,param,dir);
      _Euler2DLeftEigenvectors_(uavg,L,param,dir);
      _Euler2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler2DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,dir);
        eigL = D[k*_MODEL_NVARS_+k];
        _Euler2DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,dir);
        eigR = D[k*_MODEL_NVARS_+k];
        _Euler2DEigenvalues_(uavg,D,param,dir);
        eigC = D[k*_MODEL_NVARS_+k];

        double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
