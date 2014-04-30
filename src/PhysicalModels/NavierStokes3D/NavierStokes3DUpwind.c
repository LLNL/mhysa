#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

int NavierStokes3DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar           *solver = (HyPar*)    s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->physics;
  int             done;

  int *dim  = solver->dim_local;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (uR[_MODEL_NVARS_*p+4] - uL[_MODEL_NVARS_*p+4]);

      _NavierStokes3DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _NavierStokes3DEigenvalues_(uavg,D,param,dir);
      _NavierStokes3DLeftEigenvectors_(uavg,L,param,dir);
      _NavierStokes3DRightEigenvectors_(uavg,R,param,dir);

      /* Harten's Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

      MatMult5(_MODEL_NVARS_,DL,D,L);
      MatMult5(_MODEL_NVARS_,modA,R,DL);
      MatVecMult5(_MODEL_NVARS_,udiss,modA,udiff);
      
      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5 * (fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - udiss[3];
      fI[_MODEL_NVARS_*p+4] = 0.5 * (fL[_MODEL_NVARS_*p+4]+fR[_MODEL_NVARS_*p+4]) - udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int NavierStokes3DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar           *solver = (HyPar*)    s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->physics;
  int             done,k;

  int *dim  = solver->dim_local;

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

      /* Roe-Fixed upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _NavierStokes3DEigenvalues_(uavg,D,param,dir);
      _NavierStokes3DLeftEigenvectors_(uavg,L,param,dir);
      _NavierStokes3DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[_MODEL_NVARS_],eigC[_MODEL_NVARS_],eigR[_MODEL_NVARS_];
      _NavierStokes3DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,dir); 
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _NavierStokes3DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,dir); 
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _NavierStokes3DEigenvalues_(uavg,D,param,dir); 
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

int NavierStokes3DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar           *solver = (HyPar*)    s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->physics;
  int             done;

  int *dim  = solver->dim_local;

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

      /* Local Lax-Friedrich upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _NavierStokes3DEigenvalues_(uavg,D,param,dir);
      _NavierStokes3DLeftEigenvectors_(uavg,L,param,dir);
      _NavierStokes3DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[_MODEL_NVARS_],eigC[_MODEL_NVARS_],eigR[_MODEL_NVARS_];
      _NavierStokes3DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,dir); 
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _NavierStokes3DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,dir); 
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _NavierStokes3DEigenvalues_(uavg,D,param,dir); 
      eigC[0] = D[0];
      eigC[1] = D[6];
      eigC[2] = D[12];
      eigC[3] = D[18];
      eigC[4] = D[24];

      double alpha;
      alpha = max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));
      alpha = max3(absolute(eigL[4]),absolute(eigC[4]),absolute(eigR[4]));
      fc[4] = 0.5 * (fcL[4] + fcR[4] + alpha * (ucL[4]-ucR[4]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult5(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}
