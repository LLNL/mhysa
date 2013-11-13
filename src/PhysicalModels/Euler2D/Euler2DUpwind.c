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

      _Euler2DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _Euler2DEigenvalues_(uavg,D,param,dir);
      _Euler2DLeftEigenvectors_(uavg,L,param,dir);
      _Euler2DRightEigenvectors_(uavg,R,param,dir);

      for (k=0; k<nvars; k++) D[k*nvars+k] = absolute(D[k*nvars+k]);
      MatMult(nvars,DL,D,L);
      MatMult(nvars,modA,R,DL);

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

      _Euler2DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _Euler2DEigenvalues_(uavg,D,param,dir);
      _Euler2DLeftEigenvectors_ (uavg,L,param,dir);
      _Euler2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult(nvars,ucL,L,(uL+nvars*p));
      MatVecMult(nvars,ucR,L,(uR+nvars*p));
      MatVecMult(nvars,fcL,L,(fL+nvars*p));
      MatVecMult(nvars,fcR,L,(fR+nvars*p));

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        _Euler2DEigenvalues_((uL+nvars*p),D,param,dir);
        eigL = D[k*nvars+k];
        _Euler2DEigenvalues_((uR+nvars*p),D,param,dir);
        eigR = D[k*nvars+k];
        _Euler2DEigenvalues_(uavg        ,D,param,dir);
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

int Euler2DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar    *solver = (HyPar*)    s;
  Euler2D  *param  = (Euler2D*)  solver->physics;
  int      done,k;

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

      _Euler2DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _Euler2DEigenvalues_(uavg,D,param,dir);
      _Euler2DLeftEigenvectors_(uavg,L,param,dir);
      _Euler2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult(nvars,ucL,L,(uL+nvars*p));
      MatVecMult(nvars,ucR,L,(uR+nvars*p));
      MatVecMult(nvars,fcL,L,(fL+nvars*p));
      MatVecMult(nvars,fcR,L,(fR+nvars*p));

      for (k = 0; k < nvars; k++) {
        double eigL,eigC,eigR;
        _Euler2DEigenvalues_((uL+nvars*p),D,param,dir);
        eigL = D[k*nvars+k];
        _Euler2DEigenvalues_((uR+nvars*p),D,param,dir);
        eigR = D[k*nvars+k];
        _Euler2DEigenvalues_(uavg,D,param,dir);
        eigC = D[k*nvars+k];

        double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult(nvars,(fI+nvars*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
