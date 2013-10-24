#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

inline int NavierStokes3DGetFlowVar (double*,double*,double*,double*,double*,double*,double*,void*);

inline int NavierStokes3DEigenvalues(double *u,double **D,void *p,int dir)
{
  NavierStokes3D  *param  = (NavierStokes3D*)  p;
  int             ierr    = 0;
  double          gamma   = param->gamma;
  double          rho,vx,vy,vz,e,P,c,vn;

  ierr = NavierStokes3DGetFlowVar(u,&rho,&vx,&vy,&vz,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  if      (dir == _XDIR_) vn = vx;
  else if (dir == _YDIR_) vn = vy;
  else if (dir == _ZDIR_) vn = vz;
  else                    vn = 0.0;

  D[0][0] = vn;     D[0][1] = 0;    D[0][2] = 0;      D[0][3] = 0;    D[0][4] = 0;
  D[1][0] = 0;      D[1][1] = vn-c; D[1][2] = 0;      D[1][3] = 0;    D[1][4] = 0;
  D[2][0] = 0;      D[2][1] = 0;    D[2][2] = vn;     D[2][3] = 0;    D[2][4] = 0;
  D[3][0] = 0;      D[3][1] = 0;    D[3][2] = 0;      D[3][3] = vn+c; D[3][4] = 0;
  D[4][0] = 0;      D[4][1] = 0;    D[4][2] = 0;      D[4][3] = 0;    D[4][4] = vn;

  return(0);
}

inline int NavierStokes3DLeftEigenvectors(double *u,double **L,void *p,int dir)
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma, K=gamma-1.0;
  double  rho,vx,vy,vz,e,P,c,csq;
  double  nx,ny,nz,lx,ly,lz,mx,my,mz;
  double  qsq,qn,qm,ql;

  ierr = NavierStokes3DGetFlowVar(u,&rho,&vx,&vy,&vz,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho); csq = c*c;
  
  if (dir == _XDIR_) {
    nx = 1.0; ny = 0.0; nz = 0.0;
    lx = 0.0; ly = 1.0; lz = 0.0;
    mx = 0.0; my = 0.0; mz = 1.0;
  } else if (dir == _YDIR_) {
    nx = 0.0; ny = 1.0; nz = 0.0;
    lx = 0.0; ly = 0.0; lz = 1.0;
    mx = 1.0; my = 0.0; mz = 0.0;
  } else if (dir == _ZDIR_) {
    nx = 0.0; ny = 0.0; nz = 1.0;
    lx = 1.0; ly = 0.0; lz = 0.0;
    mx = 0.0; my = 1.0; mz = 0.0;
  } else {
    nx = 0.0; ny = 0.0; nz = 0.0;
    lx = 0.0; ly = 0.0; lz = 0.0;
    mx = 0.0; my = 0.0; mz = 0.0;
  }
  qn  = vx*nx + vy*ny + vz*nz;
  ql  = vx*lx + vy*ly + vz*lz;
  qm  = vx*mx + vy*my + vz*mz;
  qsq = vx*vx + vy*vy + vz*vz;

  L[0][0] = (K*qsq)/(4*csq) + qn/(2*c);
  L[0][1] = - (K/(2*csq)*vx + nx/(2*c));
  L[0][2] = - (K/(2*csq)*vy + ny/(2*c));
  L[0][3] = - (K/(2*csq)*vz + nz/(2*c));
  L[0][4] = K/(2*csq);
  L[1][0] = 1.0 - (K*qsq)/(2*csq);
  L[1][1] = K*vx/csq;
  L[1][2] = K*vy/csq;
  L[1][3] = K*vz/csq;
  L[1][4] = - K/csq;
  L[2][0] = (K*qsq)/(4*csq) - qn/(2*c);
  L[2][1] = - (K/(2*csq)*vx - nx/(2*c));
  L[2][2] = - (K/(2*csq)*vy - ny/(2*c));
  L[2][3] = - (K/(2*csq)*vz - nz/(2*c));
  L[2][4] = K/(2*csq);
  L[3][0] = -ql;
  L[3][1] = lx;
  L[3][2] = ly;
  L[3][3] = lz;
  L[3][4] = 0;
  L[4][0] = -qm;
  L[4][1] = mx;
  L[4][2] = my;
  L[4][3] = mz;
  L[4][4] = 0;

  return(0);
}

inline int NavierStokes3DRightEigenvectors(double *u,double **R,void *p,int dir)
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma;
  double  rho,vx,vy,vz,e,P,c;
  double  nx,ny,nz,lx,ly,lz,mx,my,mz;
  double  qsq,qn,qm,ql,H;

  ierr = NavierStokes3DGetFlowVar(u,&rho,&vx,&vy,&vz,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);
  H    = (e+P)/rho;

  if (dir == _XDIR_) {
    nx = 1.0; ny = 0.0; nz = 0.0;
    lx = 0.0; ly = 1.0; lz = 0.0;
    mx = 0.0; my = 0.0; mz = 1.0;
  } else if (dir == _YDIR_) {
    nx = 0.0; ny = 1.0; nz = 0.0;
    lx = 0.0; ly = 0.0; lz = 1.0;
    mx = 1.0; my = 0.0; mz = 0.0;
  } else if (dir == _ZDIR_) {
    nx = 0.0; ny = 0.0; nz = 1.0;
    lx = 1.0; ly = 0.0; lz = 0.0;
    mx = 0.0; my = 1.0; mz = 0.0;
  } else {
    nx = 0.0; ny = 0.0; nz = 0.0;
    lx = 0.0; ly = 0.0; lz = 0.0;
    mx = 0.0; my = 0.0; mz = 0.0;
  }
  qn  = vx*nx + vy*ny + vz*nz;
  ql  = vx*lx + vy*ly + vz*lz;
  qm  = vx*mx + vy*my + vz*mz;
  qsq = vx*vx + vy*vy + vz*vz;

  R[0][0] = 1.0;
  R[0][1] = 1.0;
  R[0][2] = 1.0;
  R[0][3] = 0.0;
  R[0][4] = 0.0;
  R[1][0] = vx - c*nx;;
  R[1][1] = vx;
  R[1][2] = vx + c*nx;
  R[1][3] = lx;
  R[1][4] = mx;
  R[2][0] = vy - c*ny;
  R[2][1] = vy;
  R[2][2] = vy + c*ny;
  R[2][3] = ly;
  R[2][4] = my;
  R[3][0] = vz - c*nz;
  R[3][1] = vz;
  R[3][2] = vz + c*nz;
  R[3][3] = lz;
  R[3][4] = mz;
  R[4][0] = H - qn*c;
  R[4][1] = qsq / 2.0;
  R[4][2] = H + qn*c;
  R[4][3] = ql;
  R[4][4] = qm;

  return(0);
}
