#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DGetFlowVar (double*,double*,double*,double*,double*,void*);

int Euler1DEigenvalues(double *u,double *D,void *p,int dir)
{
  Euler1D *param  = (Euler1D*)  p;
  double  gamma   = param->gamma;
  double  rho,v,e,P,c;
  _DECLARE_IERR_;

  IERR Euler1DGetFlowVar(u,&rho,&v,&e,&P,param); CHECKERR(ierr);
  c = sqrt(gamma*P/rho);

  D[0*_MODEL_NVARS_+0] = (v-c);  D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;
  D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = v;    D[1*_MODEL_NVARS_+2] = 0;
  D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = (v+c);

  return(0);
}

int Euler1DLeftEigenvectors(double *u,double *L,void *p,int dir)
{
  Euler1D *param  = (Euler1D*)  p;
  double  gamma   = param->gamma;
  double  rho,v,e,P,c;
  _DECLARE_IERR_;

  IERR Euler1DGetFlowVar(u,&rho,&v,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  L[0*_MODEL_NVARS_+0] = ((gamma - 1)/(rho*c)) * (-(v*v)/2 - c*v/(gamma-1));
  L[0*_MODEL_NVARS_+1] = ((gamma - 1)/(rho*c)) * (v + c/(gamma-1));
  L[0*_MODEL_NVARS_+2] = ((gamma - 1)/(rho*c)) * (-1);
  L[1*_MODEL_NVARS_+0] = ((gamma - 1)/(rho*c)) * (rho*(-(v*v)/2+c*c/(gamma-1))/c);
  L[1*_MODEL_NVARS_+1] = ((gamma - 1)/(rho*c)) * (rho*v/c);
  L[1*_MODEL_NVARS_+2] = ((gamma - 1)/(rho*c)) * (-rho/c);
  L[2*_MODEL_NVARS_+0] = ((gamma - 1)/(rho*c)) * ((v*v)/2 - c*v/(gamma-1));
  L[2*_MODEL_NVARS_+1] = ((gamma - 1)/(rho*c)) * (-v + c/(gamma-1));
  L[2*_MODEL_NVARS_+2] = ((gamma - 1)/(rho*c)) * (1);

  return(0);
}

int Euler1DRightEigenvectors(double *u,double *R,void *p,int dir)
{
  Euler1D *param  = (Euler1D*)  p;
  double  gamma   = param->gamma;
  double  rho,v,e,P,c;
  _DECLARE_IERR_;

  IERR Euler1DGetFlowVar(u,&rho,&v,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  R[0*_MODEL_NVARS_+0] = - rho/(2*c);  R[1*_MODEL_NVARS_+0] = -rho*(v-c)/(2*c); R[2*_MODEL_NVARS_+0] = -rho*((v*v)/2+(c*c)/(gamma-1)-c*v)/(2*c);
  R[0*_MODEL_NVARS_+1] = 1;            R[1*_MODEL_NVARS_+1] = v;                R[2*_MODEL_NVARS_+1] = v*v / 2;
  R[0*_MODEL_NVARS_+2] = rho/(2*c);    R[1*_MODEL_NVARS_+2] = rho*(v+c)/(2*c);  R[2*_MODEL_NVARS_+2] = rho*((v*v)/2+(c*c)/(gamma-1)+c*v)/(2*c);

  return(0);
}
