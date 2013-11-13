#include <math.h>
#include <basic.h>
#include <physicalmodels/euler2d.h>

int Euler2DGetFlowVar(double *u,double *rho,double *vx,double *vy,
                             double *e,double *P,void *p)
{
  Euler2D *param  = (Euler2D*) p;
  double  gamma   = param->gamma, vsq;

  *rho = u[0];
  *vx  = u[1] / *rho;
  *vy  = u[2] / *rho;
  *e   = u[3];
  vsq  = (*vx)*(*vx) + (*vy)*(*vy);
  *P   = ((*e) - 0.5*(*rho)*vsq) * (gamma-1.0);

  return(0);
}

int Euler2DSetFlux(double *f,double rho,double vx,double vy,
                          double e,double P,void *p,int dir)
{
  if (dir == _XDIR_) {
    f[0] = rho * vx;
    f[1] = rho * vx * vx + P;
    f[2] = rho * vx * vy;
    f[3] = (e + P) * vx;
  } else if (dir == _YDIR_) {
    f[0] = rho * vy;
    f[1] = rho * vy * vx;
    f[2] = rho * vy * vy + P;
    f[3] = (e + P) * vy;
  }
  return(0);
}

int Euler2DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  Euler2D *param  = (Euler2D*) p;
  double  rho ,vx, vy, e ,P ,H ,csq, vsq;
  double  rhoL,vxL,vyL,eL,PL,HL,cLsq;
  double  rhoR,vxR,vyR,eR,PR,HR,cRsq;
  double  gamma = param->gamma;
  _DECLARE_IERR_;

  IERR Euler2DGetFlowVar(uL,&rhoL,&vxL,&vyL,&eL,&PL,param); CHECKERR(ierr);
  cLsq = gamma * PL/rhoL;
  HL = 0.5*(vxL*vxL+vyL*vyL) + cLsq / (gamma-1.0);

  IERR Euler2DGetFlowVar(uR,&rhoR,&vxR,&vyR,&eR,&PR,param); CHECKERR(ierr);
  cRsq = gamma * PR/rhoR;
  HR = 0.5*(vxR*vxR+vyR*vyR) + cRsq / (gamma-1.0);

  double tL = sqrt(rhoL);
  double tR = sqrt(rhoR);
  
  rho = tL * tR;
  vx  = (tL*vxL + tR*vxR) / (tL + tR);
  vy  = (tL*vyL + tR*vyR) / (tL + tR);
  H   = (tL*HL + tR*HR) / (tL + tR);
  vsq = vx*vx + vy*vy;
  csq = (gamma-1.0) * (H-0.5*vsq);
  P   = csq * rho / gamma;
  e   = P/(gamma-1.0) + 0.5*rho*vsq;

  uavg[0] = rho;
  uavg[1] = rho*vx;
  uavg[2] = rho*vy;
  uavg[3] = e;

  return(0);
}
