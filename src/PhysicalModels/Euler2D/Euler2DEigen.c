#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

/*

  Rohde, A., "Eigenvalues and Eigenvectors of the Euler Equations
  in General Geometries", AIAA 2001-2609

*/

inline int Euler2DGetFlowVar (double*,double*,double*,double*,double*,double*,void*);

inline int Euler2DEigenvalues(double *u,double **D,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  int     ierr    = 0;
  double  gamma   = param->gamma;
  double  rho,vx,vy,e,P,c,vn;

  ierr = Euler2DGetFlowVar(u,&rho,&vx,&vy,&e,&P,param); CHECKERR(ierr);
  c    = sqrt(gamma*P/rho);

  if      (dir == _XDIR_) vn = vx;
  else if (dir == _YDIR_) vn = vy;
  else                    vn = 0.0;

  D[0][0] = vn-c;   D[0][1] = 0;    D[0][2] = 0;      D[0][3] = 0;
  D[1][0] = 0;      D[1][1] = vn;   D[1][2] = 0;      D[1][3] = 0;
  D[2][0] = 0;      D[2][1] = 0;    D[2][2] = vn+c;   D[2][3] = 0;
  D[3][0] = 0;      D[3][1] = 0;    D[3][2] = 0;      D[3][3] = vn;

  return(0);
}

inline int Euler2DLeftEigenvectors(double *u,double **L,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  int     ierr  = 0;
  double  ga    = param->gamma, ga_minus_one=ga-1.0;
  double  rho,vx,vy,e,P,a,un,ek;
  double  nx = 0,ny = 0;

  ierr = Euler2DGetFlowVar(u,&rho,&vx,&vy,&e,&P,param); CHECKERR(ierr);
	ek = 0.5 * (vx*vx + vy*vy);
	a = sqrt(ga * P / rho);

  if (dir == _XDIR_) { 
    un = vx; 
    nx = 1.0; 

		L[0][0] = (ga_minus_one*ek + a*un) / (2*a*a);
		L[0][1] = ((-ga_minus_one)*vx - a*nx) / (2*a*a);
		L[0][2] = ((-ga_minus_one)*vy - a*ny) / (2*a*a);
		L[0][3] = ga_minus_one / (2*a*a);

		L[1][0] = (a*a - ga_minus_one*ek) / (a*a);
		L[1][1] = (ga_minus_one*vx) / (a*a);
		L[1][2] = (ga_minus_one*vy) / (a*a);
		L[1][3] = (-ga_minus_one) / (a*a);

		L[2][0] = (ga_minus_one*ek - a*un) / (2*a*a);
		L[2][1] = ((-ga_minus_one)*vx + a*nx) / (2*a*a);
		L[2][2] = ((-ga_minus_one)*vy + a*ny) / (2*a*a);
		L[2][3] = ga_minus_one / (2*a*a);

		L[3][0] = (vy - un*ny) / nx;
		L[3][1] = ny;
		L[3][2] = (ny*ny - 1.0) / nx;
		L[3][3] = 0.0;

  } else if (dir == _YDIR_) { 
    un = vy; 
    ny = 1.0; 

		L[0][0] = (ga_minus_one*ek+a*un) / (2*a*a);
		L[0][1] = ((1.0-ga)*vx - a*nx) / (2*a*a);
		L[0][2] = ((1.0-ga)*vy - a*ny) / (2*a*a);
		L[0][3] = ga_minus_one / (2*a*a);
 
		L[1][0] = (a*a-ga_minus_one*ek) / (a*a);
		L[1][1] = ga_minus_one*vx / (a*a);
		L[1][2] = ga_minus_one*vy / (a*a);
		L[1][3] = (1.0 - ga) / (a*a);
 
		L[2][0] = (ga_minus_one*ek-a*un) / (2*a*a);
		L[2][1] = ((1.0-ga)*vx + a*nx) / (2*a*a);
		L[2][2] = ((1.0-ga)*vy + a*ny) / (2*a*a);
		L[2][3] = ga_minus_one / (2*a*a);

		L[3][0] = (un*nx-vx) / ny;
		L[3][1] = (1.0 - nx*nx) / ny;
		L[3][2] = - nx;
		L[3][3] = 0;

  } else {
    fprintf(stderr,"Error in Euler2DLeftEigenvectors(): invalid dir!\n");
    return(1);
  }

  return(0);
}

inline int Euler2DRightEigenvectors(double *u,double **R,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  int     ierr = 0;
  double  ga   = param->gamma, ga_minus_one = ga-1.0;
  double  rho,vx,vy,e,P,un,ek,a,h0;
  double  nx = 0,ny = 0;

  ierr = Euler2DGetFlowVar(u,&rho,&vx,&vy,&e,&P,param); CHECKERR(ierr);
	ek   = 0.5 * (vx*vx + vy*vy);
	a    = sqrt(ga * P / rho);
  h0   = a*a / ga_minus_one + ek;

	if (dir == _XDIR_) {
  	un = vx;
    nx = 1.0;

		R[0][0] = 1.0;
		R[1][0] = vx - a*nx;
		R[2][0] = vy - a*ny;
		R[3][0] = h0 - a*un;

		R[0][1] = 1.0;
		R[1][1] = vx;
		R[2][1] = vy;
		R[3][1] = ek;

		R[0][2] = 1.0;
		R[1][2] = vx + a*nx;
		R[2][2] = vy + a*ny;
		R[3][2] = h0 + a*un;

		R[0][3] = 0.0;
		R[1][3] = ny;
		R[2][3] = -nx;
		R[3][3] = vx*ny - vy*nx;

	} else if (dir == _YDIR_) {
    un = vy;
    ny = 1.0;

		R[0][0] = 1.0;
		R[1][0] = vx - a*nx;
		R[2][0] = vy - a*ny;
		R[3][0] = h0 - a*un;

		R[0][1] = 1.0;
		R[1][1] = vx;
		R[2][1] = vy;
		R[3][1] = ek;

		R[0][2] = 1.0;
		R[1][2] = vx + a*nx;
		R[2][2] = vy + a*ny;
		R[3][2] = h0 + a*un;

		R[0][3] = 0;
		R[1][3] = ny;
		R[2][3] = -nx;
		R[3][3] = vx*ny-vy*nx;

  } else {
    fprintf(stderr,"Error in Euler2DRightEigenvectors(): invalid dir!\n");
    return(1);
  }

  return(0);
}
