#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DSourceUpwindRF(double *fI,double *fL,double *fR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_],
             ugL[_MODEL_NVARS_], ugR[_MODEL_NVARS_];
      double kappa = max(param->grav_field[pL],param->grav_field[pR]);

      /* Roe-Fixed upwinding scheme */

      _Euler1DRoeAverage_       (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _Euler1DEigenvalues_      (uavg,D,param,0);
      _Euler1DLeftEigenvectors_ (uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* gravity-corrected state vector */
      _ArrayCopy1D_((u+_MODEL_NVARS_*pL),ugL,_MODEL_NVARS_); _ArrayScale1D_(ugL,(1.0/param->grav_field[pL]),_MODEL_NVARS_);
      _ArrayCopy1D_((u+_MODEL_NVARS_*pR),ugR,_MODEL_NVARS_); _ArrayScale1D_(ugR,(1.0/param->grav_field[pR]),_MODEL_NVARS_);

      /* calculate characteristic fluxes and variables */
      MatVecMult3(_MODEL_NVARS_,ucL,L,ugL);
      MatVecMult3(_MODEL_NVARS_,ucR,L,ugR);
      MatVecMult3(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((u+_MODEL_NVARS_*pL),D,param,0);
        eigL = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_((u+_MODEL_NVARS_*pR),D,param,0);
        eigR = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = D[k*_MODEL_NVARS_+k];

        if ((eigL > 0) && (eigC > 0) && (eigR > 0))       fc[k] = fcL[k];
        else if ((eigL < 0) && (eigC < 0) && (eigR < 0))  fc[k] = fcR[k];
        else {
          double alpha = kappa * max3(absolute(eigL),absolute(eigC),absolute(eigR));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult3(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler1DSourceUpwindLLF(double *fI,double *fL,double *fR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_],
             ugL[_MODEL_NVARS_], ugR[_MODEL_NVARS_];
      double kappa = max(param->grav_field[pL],param->grav_field[pR]);

      /* Local Lax-Friedrich upwinding scheme */

      _Euler1DRoeAverage_       (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _Euler1DEigenvalues_      (uavg,D,param,0);
      _Euler1DLeftEigenvectors_ (uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* gravity-corrected state vector */
      _ArrayCopy1D_((u+_MODEL_NVARS_*pL),ugL,_MODEL_NVARS_); _ArrayScale1D_(ugL,(1.0/param->grav_field[pL]),_MODEL_NVARS_);
      _ArrayCopy1D_((u+_MODEL_NVARS_*pR),ugR,_MODEL_NVARS_); _ArrayScale1D_(ugR,(1.0/param->grav_field[pR]),_MODEL_NVARS_);

      /* calculate characteristic fluxes and variables */
      MatVecMult3(_MODEL_NVARS_,ucL,L,ugL);
      MatVecMult3(_MODEL_NVARS_,ucR,L,ugR);
      MatVecMult3(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((u+_MODEL_NVARS_*pL),D,param,0);
        eigL = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_((u+_MODEL_NVARS_*pR),D,param,0);
        eigR = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = D[k*_MODEL_NVARS_+k];

        double alpha = kappa * max3(absolute(eigL),absolute(eigC),absolute(eigR));
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult3(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
