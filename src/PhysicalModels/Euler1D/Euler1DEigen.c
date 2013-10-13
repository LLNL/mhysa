#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

inline int Euler1DGetFlowVar (double*,double*,double*,double*,double*,void*);

inline int Euler1DEigenvalues(double *u,double **D,void *p)
{
  Euler1D *param  = (Euler1D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma;
  double  rho,v,e,P,c;

  ierr = Euler1DGetFlowVar(u,&rho,&v,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  D[0][0] = (v-c);  D[0][1] = 0;    D[0][2] = 0;
  D[1][0] = 0;      D[1][1] = v;    D[1][2] = 0;
  D[2][0] = 0;      D[2][1] = 0;    D[2][2] = (v+c);

  return(0);
}

inline int Euler1DLeftEigenvectors(double *u,double **L,void *p)
{
  Euler1D *param  = (Euler1D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma;
  double  rho,v,e,P,c;

  ierr = Euler1DGetFlowVar(u,&rho,&v,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  L[0][0] = ((gamma - 1)/(rho*c)) * (-(v*v)/2 - c*v/(gamma-1));
  L[0][1] = ((gamma - 1)/(rho*c)) * (v + c/(gamma-1));
  L[0][2] = ((gamma - 1)/(rho*c)) * (-1);
  L[1][0] = ((gamma - 1)/(rho*c)) * (rho*(-(v*v)/2+c*c/(gamma-1))/c);
  L[1][1] = ((gamma - 1)/(rho*c)) * (rho*v/c);
  L[1][2] = ((gamma - 1)/(rho*c)) * (-rho/c);
  L[2][0] = ((gamma - 1)/(rho*c)) * ((v*v)/2 - c*v/(gamma-1));
  L[2][1] = ((gamma - 1)/(rho*c)) * (-v + c/(gamma-1));
  L[2][2] = ((gamma - 1)/(rho*c)) * (1);

  return(0);
}

inline int Euler1DRightEigenvectors(double *u,double **R,void *p)
{
  Euler1D *param  = (Euler1D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma;
  double  rho,v,e,P,c;

  ierr = Euler1DGetFlowVar(u,&rho,&v,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  R[0][0] = - rho/(2*c);  R[1][0] = -rho*(v-c)/(2*c); R[2][0] = -rho*((v*v)/2+(c*c)/(gamma-1)-c*v)/(2*c);
  R[0][1] = 1;            R[1][1] = v;                R[2][1] = v*v / 2;
  R[0][2] = rho/(2*c);    R[1][2] = rho*(v+c)/(2*c);  R[2][2] = rho*((v*v)/2+(c*c)/(gamma-1)+c*v)/(2*c);

  return(0);
}
