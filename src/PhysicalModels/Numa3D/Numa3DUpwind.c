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

  int *dim  = solver->dim_local;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_];

      /* Rusanov's upwinding scheme */

      uavg[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] + uL[_MODEL_NVARS_*p+0]);
      uavg[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] + uL[_MODEL_NVARS_*p+1]);
      uavg[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] + uL[_MODEL_NVARS_*p+2]);
      uavg[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] + uL[_MODEL_NVARS_*p+3]);
      uavg[4] = 0.5 * (uR[_MODEL_NVARS_*p+4] + uL[_MODEL_NVARS_*p+4]);

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (uR[_MODEL_NVARS_*p+4] - uL[_MODEL_NVARS_*p+4]);

      double drho,vel[3],dT,dP,rho0,P0,T0,c;
      rho0 = (dir==_ZDIR_ ? 0.5*(param->rho0[index_inter[_ZDIR_]]+param->rho0[index_inter[_ZDIR_]-1]) 
                          : param->rho0[index_inter[_ZDIR_]] );
      T0   = (dir==_ZDIR_ ? 0.5*(param->T0[index_inter[_ZDIR_]]+param->T0[index_inter[_ZDIR_]-1]) 
                          : param->T0[index_inter[_ZDIR_]] );
      P0   = (dir==_ZDIR_ ? 0.5*(param->P0[index_inter[_ZDIR_]]+param->P0[index_inter[_ZDIR_]-1]) 
                          : param->P0[index_inter[_ZDIR_]] );

      _Numa3DGetFlowVars_         (uavg,drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputePressure_     (param,T0,dT,P0,dP);
      _Numa3DComputeSpeedofSound_ (param->gamma,P0,dP,rho0,drho,c);

      double alpha  =   c + absolute(vel[dir]);

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
