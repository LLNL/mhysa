#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

int Numa3DRusanov(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar   *solver = (HyPar*)  s;
  Numa3D  *param  = (Numa3D*) solver->physics;
  int      done;

  int *dim   = solver->dim_local;
  int ghosts = solver->ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_],drho,vel[3],dT,rho0,T0,EP,c;

      /* Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (uR[_MODEL_NVARS_*p+4] - uL[_MODEL_NVARS_*p+4]);

      /* left of the interface */
      rho0 = param->rho0  [(dir==_ZDIR_?index_inter[_ZDIR_]-1:index_inter[_ZDIR_])+ghosts];
      T0   = param->T0    [(dir==_ZDIR_?index_inter[_ZDIR_]-1:index_inter[_ZDIR_])+ghosts];
      EP   = param->ExnerP[(dir==_ZDIR_?index_inter[_ZDIR_]-1:index_inter[_ZDIR_])+ghosts];
      _Numa3DGetFlowVars_         ((uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeSpeedofSound_ (param->gamma,param->R,T0,dT,rho0,drho,EP,c);
      double alphaL = c + absolute(vel[dir]);

      /* right of the interface */
      rho0 = param->rho0  [index_inter[_ZDIR_]+ghosts];
      T0   = param->T0    [index_inter[_ZDIR_]+ghosts];
      EP   = param->ExnerP[index_inter[_ZDIR_]+ghosts];
      _Numa3DGetFlowVars_         ((uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeSpeedofSound_ (param->gamma,param->R,T0,dT,rho0,drho,EP,c);
      double alphaR = c + absolute(vel[dir]);

      double alpha = max(alphaL,alphaR);
      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - alpha*udiff[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - alpha*udiff[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - alpha*udiff[2];
      fI[_MODEL_NVARS_*p+3] = 0.5 * (fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - alpha*udiff[3];
      fI[_MODEL_NVARS_*p+4] = 0.5 * (fL[_MODEL_NVARS_*p+4]+fR[_MODEL_NVARS_*p+4]) - alpha*udiff[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Numa3DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar   *solver = (HyPar*)    s;
  Numa3D  *param  = (Numa3D*)  solver->physics;
  int     done,k;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      double rho0L, T0L, EPL;
      rho0L = param->rho0  [(dir==_ZDIR_?index_inter[_ZDIR_]-1:index_inter[_ZDIR_])+ghosts];
      T0L   = param->T0    [(dir==_ZDIR_?index_inter[_ZDIR_]-1:index_inter[_ZDIR_])+ghosts];
      EPL   = param->ExnerP[(dir==_ZDIR_?index_inter[_ZDIR_]-1:index_inter[_ZDIR_])+ghosts];

      double rho0R, T0R, EPR;
      rho0R = param->rho0  [index_inter[_ZDIR_]+ghosts];
      T0R   = param->T0    [index_inter[_ZDIR_]+ghosts];
      EPR   = param->ExnerP[index_inter[_ZDIR_]+ghosts];

      double rho0, T0, EP;
      rho0  = 0.5 * (rho0L + rho0R);
      T0    = 0.5 * (T0L   + T0R  );
      EP    = 0.5 * (EPL   + EPR  );

      /* Roe-Fixed upwinding scheme */

      _Numa3DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param,
                         rho0L,rho0R,rho0,T0L,T0R,T0,EPL,EPR,EP);

      _Numa3DLeftEigenvectors_ (uavg,L,param,dir,rho0,T0,EP);
      _Numa3DRightEigenvectors_(uavg,R,param,dir,rho0,T0,EP);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[_MODEL_NVARS_],eigC[_MODEL_NVARS_],eigR[_MODEL_NVARS_];
      _Numa3DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,dir,rho0,T0,EP); 
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _Numa3DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,dir,rho0,T0,EP); 
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _Numa3DEigenvalues_(uavg,D,param,dir,rho0,T0,EP); 
      eigC[0] = D[0];
      eigC[1] = D[6];
      eigC[2] = D[12];
      eigC[3] = D[18];
      eigC[4] = D[24];

      for (k = 0; k < _MODEL_NVARS_; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult5(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

