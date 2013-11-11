#include <math.h>
#include <basic.h>
#include <physicalmodels/euler1d.h>

int Euler1DGetFlowVar(double *u,double *rho,double *v,double *e,double *P,void *p)
{
  Euler1D   *param  = (Euler1D*) p;
  double    gamma   = param->gamma;

  *rho = u[0];
  *v   = u[1] / *rho;
  *e   = u[2];
  *P   = ((*e) - 0.5*(*rho)*(*v)*(*v)) * (gamma-1.0);

  return(0);
}

int Euler1DSetFlux(double *f,double rho,double v,double e,double P,void *p)
{
  f[0] = rho * v;
  f[1] = rho * v * v + P;
  f[2] = (e + P) * v;

  return(0);
}

int Euler1DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  Euler1D *param  = (Euler1D*) p;
  int     ierr = 0;
  double  rho ,v ,e ,P ,H ,csq;
  double  rhoL,vL,eL,PL,HL,cLsq;
  double  rhoR,vR,eR,PR,HR,cRsq;
  double  gamma = param->gamma;

  ierr = Euler1DGetFlowVar(uL,&rhoL,&vL,&eL,&PL,param); CHECKERR(ierr);
  cLsq = gamma * PL/rhoL;
  HL = 0.5*vL*vL + cLsq / (gamma-1.0);

  ierr = Euler1DGetFlowVar(uR,&rhoR,&vR,&eR,&PR,param); CHECKERR(ierr);
  cRsq = gamma * PR/rhoR;
  HR = 0.5*vR*vR + cRsq / (gamma-1.0);

  double tL = sqrt(rhoL);
  double tR = sqrt(rhoR);
  
  rho = tL * tR;
  v   = (tL*vL + tR*vR) / (tL + tR);
  H   = (tL*HL + tR*HR) / (tL + tR);
  csq = (gamma-1.0) * (H-0.5*v*v);
  P   = csq * rho / gamma;
  e   = P/(gamma-1.0) + 0.5*rho*v*v;

  uavg[0] = rho;
  uavg[1] = rho*v;
  uavg[2] = e;

  return(0);
}
