#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d-cons.h>
#include <hypar.h>

int Numa2DConsRusanov(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar       *solver = (HyPar*)  s;
  Numa2DCons  *param  = (Numa2DCons*) solver->physics;
  int         done;

  int *dim   = solver->dim_local;
  int ghosts = solver->ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D2_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D2_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    _ArrayCopy1D2_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_],rho,vel[2],theta,c;

      /* Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      /* left of the interface */
      _Numa2DConsGetFlowVars_         ((uL+_MODEL_NVARS_*p),rho,vel[0],vel[1],theta);
      _Numa2DConsComputeSpeedofSound_ (param,rho,theta,c);
      double alphaL = c + absolute(vel[dir]);

      /* right of the interface */
      _Numa2DConsGetFlowVars_         ((uR+_MODEL_NVARS_*p),rho,vel[0],vel[1],theta);
      _Numa2DConsComputeSpeedofSound_ (param,rho,theta,c);
      double alphaR = c + absolute(vel[dir]);

      double alpha = max(alphaL,alphaR);
      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - alpha*udiff[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - alpha*udiff[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - alpha*udiff[2];
      fI[_MODEL_NVARS_*p+3] = 0.5 * (fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - alpha*udiff[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

