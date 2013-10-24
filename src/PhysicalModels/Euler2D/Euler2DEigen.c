#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

inline int Euler2DGetFlowVar (double*,double*,double*,double*,double*,double*,void*);

inline int Euler2DEigenvalues(double *u,double **D,void *p,int dir)
{
  Euler2D  *param  = (Euler2D*)  p;
  int      ierr    = 0;
  double   gamma   = param->gamma;
  double   rho,vx,vy,e,P,c,vn;

  ierr = Euler2DGetFlowVar(u,&rho,&vx,&vy,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  if      (dir == _XDIR_) vn = vx;
  else if (dir == _YDIR_) vn = vy;
  else                    vn = 0;

  D[0][0] = vn-c; D[0][1] = 0;      D[0][2] = 0;    D[0][3] = 0;
  D[1][0] = 0;    D[1][1] = vn;     D[1][2] = 0;    D[1][3] = 0;
  D[2][0] = 0;    D[2][1] = 0;      D[2][2] = vn+c; D[2][3] = 0;
  D[3][0] = 0;    D[3][1] = 0;      D[3][2] = 0;    D[3][3] = vn;

  return(0);
}

inline int Euler2DLeftEigenvectors(double *u,double **L,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma, K=gamma-1.0;
  double  rho,vx,vy,e,P,c,csq;
  double  nx,ny,lx,ly;
  double  qsq,qn,ql;

  ierr = Euler2DGetFlowVar(u,&rho,&vx,&vy,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho); csq = c*c;
  
  if (dir == _XDIR_) {
    nx = 1.0; ny = 0.0;
    lx = 0.0; ly = 1.0;
  } else if (dir == _YDIR_) {
    nx = 0.0; ny = 1.0;
    lx = 1.0; ly = 0.0;
  } else {
    nx = 0.0; ny = 0.0;
    lx = 0.0; ly = 0.0;
  }
  qn  = vx*nx + vy*ny;
  ql  = vx*lx + vy*ly;
  qsq = vx*vx + vy*vy;

  L[0][0] = (K*qsq)/(4*csq) + qn/(2*c);
  L[0][1] = - (K/(2*csq)*vx + nx/(2*c));
  L[0][2] = - (K/(2*csq)*vy + ny/(2*c));
  L[0][3] = K/(2*csq);
  L[1][0] = 1.0 - (K*qsq)/(2*csq);
  L[1][1] = K*vx/csq;
  L[1][2] = K*vy/csq;
  L[1][3] = - K/csq;
  L[2][0] = (K*qsq)/(4*csq) - qn/(2*c);
  L[2][1] = - (K/(2*csq)*vx - nx/(2*c));
  L[2][2] = - (K/(2*csq)*vy - ny/(2*c));
  L[2][3] = K/(2*csq);
  L[3][0] = -ql;
  L[3][1] = lx;
  L[3][2] = ly;
  L[3][3] = 0;

  return(0);
}

inline int Euler2DRightEigenvectors(double *u,double **R,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma;
  double  rho,vx,vy,e,P,c;
  double  nx,ny,lx,ly;
  double  qsq,qn,ql,H;

  ierr = Euler2DGetFlowVar(u,&rho,&vx,&vy,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);
  H    = (e+P)/rho;

  if (dir == _XDIR_) {
    nx = 1.0; ny = 0.0;
    lx = 0.0; ly = 1.0;
  } else if (dir == _YDIR_) {
    nx = 0.0; ny = 1.0;
    lx = 1.0; ly = 0.0;
  } else {
    nx = 0.0; ny = 0.0;
    lx = 0.0; ly = 0.0;
  }
  qn  = vx*nx + vy*ny;
  ql  = vx*lx + vy*ly;
  qsq = vx*vx + vy*vy;

  R[0][0] = 1.0;
  R[0][1] = 1.0;
  R[0][2] = 1.0;
  R[0][3] = 0.0;
  R[1][0] = vx - c*nx;;
  R[1][1] = vx;
  R[1][2] = vx + c*nx;
  R[1][3] = lx;
  R[2][0] = vy - c*ny;
  R[2][1] = vy;
  R[2][2] = vy + c*ny;
  R[2][3] = ly;
  R[3][0] = H - qn*c;
  R[3][1] = qsq / 2.0;
  R[3][2] = H + qn*c;
  R[3][3] = ql;

  return(0);
}
