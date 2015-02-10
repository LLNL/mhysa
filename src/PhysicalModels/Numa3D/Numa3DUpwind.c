#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

int Numa3DRusanovFlux(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
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
      double udiff[_MODEL_NVARS_],drho,vel[3],dT,rho0,P0,T0,EP,c;
      double zcoordL, zcoordR;

      if (dir == _ZDIR_) {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]-1),dim,ghosts,solver->x,zcoordL);
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->x,zcoordR);
      } else {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->x,zcoordL);
        zcoordR = zcoordL;
      }

      /* Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (uR[_MODEL_NVARS_*p+4] - uL[_MODEL_NVARS_*p+4]);

      /* left of the interface */
      param->StandardAtmosphere(param,zcoordL,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeSpeedofSound_ (param->gamma,param->R,T0,dT,rho0,drho,EP,c);
      double alphaL = c + absolute(vel[dir]);

      /* right of the interface */
      param->StandardAtmosphere(param,zcoordR,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeSpeedofSound_ (param->gamma,param->R,T0,dT,rho0,drho,EP,c);
      double alphaR = c + absolute(vel[dir]);

      double alpha = max(alphaL,alphaR);
      fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
      fI[_MODEL_NVARS_*p+4] = 0.5*(fL[_MODEL_NVARS_*p+4]+fR[_MODEL_NVARS_*p+4])-alpha*udiff[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Numa3DRusanovLinearFlux(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
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
      double udiff[_MODEL_NVARS_],drho,vel[3],dT,rho0,P0,T0,EP,c;
      double zcoordL, zcoordR;

      if (dir == _ZDIR_) {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]-1),dim,ghosts,solver->x,zcoordL);
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->x,zcoordR);
      } else {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->x,zcoordL);
        zcoordR = zcoordL;
      }

      /* Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (uR[_MODEL_NVARS_*p+4] - uL[_MODEL_NVARS_*p+4]);

      /* left of the interface */
      param->StandardAtmosphere   (param,zcoordL,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeLinearizedSpeedofSound_ (param->gamma,param->R,T0,rho0,EP,c);
      double alphaL = c + absolute(vel[dir]);

      /* right of the interface */
      param->StandardAtmosphere(param,zcoordR,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeLinearizedSpeedofSound_ (param->gamma,param->R,T0,rho0,EP,c);
      double alphaR = c + absolute(vel[dir]);

      double alpha = max(alphaL,alphaR);
      fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
      fI[_MODEL_NVARS_*p+4] = 0.5*(fL[_MODEL_NVARS_*p+4]+fR[_MODEL_NVARS_*p+4])-alpha*udiff[4];

      vel[0] = dT; /* useless statement to avoid compiler warning */
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}
