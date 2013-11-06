#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <physicalmodels/navierstokes3d.h>

inline int NavierStokes3DGetFlowVar(double *u,double *rho,
                                    double *vx,double *vy,double *vz,
                                    double *e,double *P,void *p)
{
  NavierStokes3D   *param  = (NavierStokes3D*) p;
  double           gamma   = param->gamma, vsq;

  *rho = u[0];
  *vx  = u[1] / *rho;
  *vy  = u[2] / *rho;
  *vz  = u[3] / *rho;
  *e   = u[4];
  vsq  = (*vx)*(*vx) + (*vy)*(*vy) + (*vz)*(*vz);
  *P   = ((*e) - 0.5*(*rho)*vsq) * (gamma-1.0);
  if (*rho < 0) {
    fprintf(stderr,"Error: Negative density encountered.\n");
    return(1);
  }
  if (*P < 0) {
    fprintf(stderr,"Error: Negative pressure encountered.\n");
    return(1);
  }

  return(0);
}

inline int NavierStokes3DSetFlux(double *f,double rho,
                                 double vx,double vy,double vz,
                                 double e,double P,void *p,int dir)
{
  if (dir == _XDIR_) {
    f[0] = rho * vx;
    f[1] = rho * vx * vx + P;
    f[2] = rho * vx * vy;
    f[3] = rho * vx * vz;
    f[4] = (e + P) * vx;
  } else if (dir == _YDIR_) {
    f[0] = rho * vy;
    f[1] = rho * vy * vx;
    f[2] = rho * vy * vy + P;
    f[3] = rho * vy * vz;
    f[4] = (e + P) * vy;
  } else if (dir == _ZDIR_) {
    f[0] = rho * vz;
    f[1] = rho * vz * vx;
    f[2] = rho * vz * vy;
    f[3] = rho * vz * vz + P;
    f[4] = (e + P) * vz;
  }
  return(0);
}

inline int NavierStokes3DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  NavierStokes3D *param  = (NavierStokes3D*) p;
  int     ierr = 0;
  double  rho ,vx, vy, vz, e ,P ,H;
  double  rhoL,vxL,vyL,vzL,eL,PL,HL,cLsq;
  double  rhoR,vxR,vyR,vzR,eR,PR,HR,cRsq;
  double  gamma = param->gamma;

  ierr = NavierStokes3DGetFlowVar(uL,&rhoL,&vxL,&vyL,&vzL,&eL,&PL,param); CHECKERR(ierr);
  cLsq = gamma * PL/rhoL;
  HL = 0.5*(vxL*vxL+vyL*vyL+vzL*vzL) + cLsq / (gamma-1.0);

  ierr = NavierStokes3DGetFlowVar(uR,&rhoR,&vxR,&vyR,&vzR,&eR,&PR,param); CHECKERR(ierr);
  cRsq = gamma * PR/rhoR;
  HR = 0.5*(vxR*vxR+vyR*vyR+vzR*vzR) + cRsq / (gamma-1.0);

  double tL = sqrt(rhoL);
  double tR = sqrt(rhoR);
  
  rho = tL * tR;
  vx  = (tL*vxL + tR*vxR) / (tL + tR);
  vy  = (tL*vyL + tR*vyR) / (tL + tR);
  vz  = (tL*vzL + tR*vzR) / (tL + tR);
  H   = (tL*HL  + tR*HR ) / (tL + tR);
  P = (H - 0.5* (vx*vx+vy*vy+vz*vz)) * (rho*(gamma-1.0))/gamma;
  e   = P/(gamma-1.0) + 0.5*rho*(vx*vx+vy*vy+vz*vz);

  uavg[0] = rho;
  uavg[1] = rho*vx;
  uavg[2] = rho*vy;
  uavg[3] = rho*vz;
  uavg[4] = e;

  return(0);
}
