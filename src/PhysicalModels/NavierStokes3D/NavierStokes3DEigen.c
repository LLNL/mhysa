#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*

  Rohde, A., "Eigenvalues and Eigenvectors of the Euler Equations
  in General Geometries", AIAA 2001-2609

*/

inline int NavierStokes3DGetFlowVar (double*,double*,double*,double*,double*,double*,double*,void*);

inline int NavierStokes3DEigenvalues(double *u,double *D,void *p,int dir)
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

  D[0*_MODEL_NVARS_+0] = vn-c;   D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0;    D[0*_MODEL_NVARS_+4] = 0;
  D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vn;   D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0;    D[1*_MODEL_NVARS_+4] = 0;
  D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vn+c;   D[2*_MODEL_NVARS_+3] = 0;    D[2*_MODEL_NVARS_+4] = 0;
  D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vn;   D[3*_MODEL_NVARS_+4] = 0;
  D[4*_MODEL_NVARS_+0] = 0;      D[4*_MODEL_NVARS_+1] = 0;    D[4*_MODEL_NVARS_+2] = 0;      D[4*_MODEL_NVARS_+3] = 0;    D[4*_MODEL_NVARS_+4] = vn;

  return(0);
}

inline int NavierStokes3DLeftEigenvectors(double *u,double *L,void *p,int dir)
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  int     ierr  = 0;
  double  ga    = param->gamma, ga_minus_one=ga-1.0;
  double  rho,vx,vy,vz,e,P,a,un,ek;
  double  nx = 0,ny = 0,nz = 0;

  ierr = NavierStokes3DGetFlowVar(u,&rho,&vx,&vy,&vz,&e,&P,param); CHECKERR(ierr);
	ek = 0.5 * (vx*vx + vy*vy + vz*vz);
	a = sqrt(ga * P / rho);

  if (dir == _XDIR_) { 
    un = vx; 
    nx = 1.0; 

		L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek + a*un) / (2*a*a);
		L[0*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx - a*nx) / (2*a*a);
		L[0*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy - a*ny) / (2*a*a);
		L[0*_MODEL_NVARS_+3] = ((-ga_minus_one)*vz - a*nz) / (2*a*a);
		L[0*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a);

		L[1*_MODEL_NVARS_+0] = (a*a - ga_minus_one*ek) / (a*a);
		L[1*_MODEL_NVARS_+1] = (ga_minus_one*vx) / (a*a);
		L[1*_MODEL_NVARS_+2] = (ga_minus_one*vy) / (a*a);
		L[1*_MODEL_NVARS_+3] = (ga_minus_one*vz) / (a*a);
		L[1*_MODEL_NVARS_+4] = (-ga_minus_one) / (a*a);

		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek - a*un) / (2*a*a);
		L[2*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx + a*nx) / (2*a*a);
		L[2*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy + a*ny) / (2*a*a);
		L[2*_MODEL_NVARS_+3] = ((-ga_minus_one)*vz + a*nz) / (2*a*a);
		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a);

		L[3*_MODEL_NVARS_+0] = (vy - un*ny) / nx;
		L[3*_MODEL_NVARS_+1] = ny;
		L[3*_MODEL_NVARS_+2] = (ny*ny - 1.0) / nx;
		L[3*_MODEL_NVARS_+3] = ny * nz / nx;
		L[3*_MODEL_NVARS_+4] = 0.0;

		L[4*_MODEL_NVARS_+0] = (un*nz - vz) / nx;
		L[4*_MODEL_NVARS_+1] = - nz;
		L[4*_MODEL_NVARS_+2] = - ny * nz / nx;
		L[4*_MODEL_NVARS_+3] = (1 - nz*nz) / nx;
		L[4*_MODEL_NVARS_+4] = 0.0;

  } else if (dir == _YDIR_) { 
    un = vy; 
    ny = 1.0; 

		L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*un) / (2*a*a);
		L[0*_MODEL_NVARS_+1] = ((1.0-ga)*vx - a*nx) / (2*a*a);
		L[0*_MODEL_NVARS_+2] = ((1.0-ga)*vy - a*ny) / (2*a*a);
		L[0*_MODEL_NVARS_+3] = ((1.0-ga)*vz - a*nz) / (2*a*a);
		L[0*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a);
 
		L[1*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a);
		L[1*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a);
		L[1*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a);
		L[1*_MODEL_NVARS_+3] = ga_minus_one*vz / (a*a);
		L[1*_MODEL_NVARS_+4] = (1.0 - ga) / (a*a);
 
		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*un) / (2*a*a);
		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx + a*nx) / (2*a*a);
		L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy + a*ny) / (2*a*a);
		L[2*_MODEL_NVARS_+3] = ((1.0-ga)*vz + a*nz) / (2*a*a);
		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a);

		L[3*_MODEL_NVARS_+0] = (un*nx-vx) / ny;
		L[3*_MODEL_NVARS_+1] = (1.0 - nx*nx) / ny;
		L[3*_MODEL_NVARS_+2] = - nx;
		L[3*_MODEL_NVARS_+3] = -nx*nz / ny;
		L[3*_MODEL_NVARS_+4] = 0;

		L[4*_MODEL_NVARS_+0] = (vz - un*nz) / ny;
		L[4*_MODEL_NVARS_+1] = nx*nz / ny;
		L[4*_MODEL_NVARS_+2] = nz;
		L[4*_MODEL_NVARS_+3] = (nz*nz - 1.0) / ny;
		L[4*_MODEL_NVARS_+4] = 0;
 
  } else if (dir == _ZDIR_) { 
    un = vz; 
    nz = 1.0; 

		L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*un) / (2*a*a);
		L[0*_MODEL_NVARS_+1] = ((1.0-ga)*vx-a*nx) / (2*a*a);
		L[0*_MODEL_NVARS_+2] = ((1.0-ga)*vy-a*ny) / (2*a*a);
		L[0*_MODEL_NVARS_+3] = ((1.0-ga)*vz-a*nz) / (2*a*a);
		L[0*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a);

		L[1*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a);
		L[1*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a);
		L[1*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a);
		L[1*_MODEL_NVARS_+3] = ga_minus_one*vz / (a*a);
		L[1*_MODEL_NVARS_+4] = (1.0-ga) / (a*a);

		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*un) / (2*a*a);
		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx+a*nx) / (2*a*a);
		L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy+a*ny) / (2*a*a);
		L[2*_MODEL_NVARS_+3] = ((1.0-ga)*vz+a*nz) / (2*a*a);
		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a);

		L[3*_MODEL_NVARS_+0] = (vx-un*nx) / nz;
		L[3*_MODEL_NVARS_+1] = (nx*nx-1.0) / nz;
		L[3*_MODEL_NVARS_+2] = nx*ny / nz;
		L[3*_MODEL_NVARS_+3] = nx;
		L[3*_MODEL_NVARS_+4] = 0;

		L[4*_MODEL_NVARS_+0] = (un*ny-vy) / nz;
		L[4*_MODEL_NVARS_+1] = -nx*ny/nz;
		L[4*_MODEL_NVARS_+2] = (1.0-ny*ny) / nz;
		L[4*_MODEL_NVARS_+3] = -ny;
		L[4*_MODEL_NVARS_+4] = 0;

  } else {
    fprintf(stderr,"Error in NavierStokes3DLeftEigenvectors(): invalid dir!\n");
    return(1);
  }

  return(0);
}

inline int NavierStokes3DRightEigenvectors(double *u,double *R,void *p,int dir)
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  int     ierr = 0;
  double  ga   = param->gamma, ga_minus_one = ga-1.0;
  double  rho,vx,vy,vz,e,P,un,ek,a,h0;
  double  nx = 0,ny = 0,nz = 0;

  ierr = NavierStokes3DGetFlowVar(u,&rho,&vx,&vy,&vz,&e,&P,param); CHECKERR(ierr);
	ek   = 0.5 * (vx*vx + vy*vy + vz*vz);
	a    = sqrt(ga * P / rho);
  h0   = a*a / ga_minus_one + ek;

	if (dir == _XDIR_) {
  	un = vx;
    nx = 1.0;

		R[0*_MODEL_NVARS_+0] = 1.0;
		R[1*_MODEL_NVARS_+0] = vx - a*nx;
		R[2*_MODEL_NVARS_+0] = vy - a*ny;
		R[3*_MODEL_NVARS_+0] = vz - a*nz;
		R[4*_MODEL_NVARS_+0] = h0 - a*un;

		R[0*_MODEL_NVARS_+1] = 1.0;
		R[1*_MODEL_NVARS_+1] = vx;
		R[2*_MODEL_NVARS_+1] = vy;
		R[3*_MODEL_NVARS_+1] = vz;
		R[4*_MODEL_NVARS_+1] = ek;

		R[0*_MODEL_NVARS_+2] = 1.0;
		R[1*_MODEL_NVARS_+2] = vx + a*nx;
		R[2*_MODEL_NVARS_+2] = vy + a*ny;
		R[3*_MODEL_NVARS_+2] = vz + a*nz;
		R[4*_MODEL_NVARS_+2] = h0 + a*un;

		R[0*_MODEL_NVARS_+3] = 0.0;
		R[1*_MODEL_NVARS_+3] = ny;
		R[2*_MODEL_NVARS_+3] = -nx;
		R[3*_MODEL_NVARS_+3] = 0.0;
		R[4*_MODEL_NVARS_+3] = vx*ny - vy*nx;

		R[0*_MODEL_NVARS_+4] = 0.0;
		R[1*_MODEL_NVARS_+4] = -nz;
		R[2*_MODEL_NVARS_+4] = 0.0;
		R[3*_MODEL_NVARS_+4] = nx;
		R[4*_MODEL_NVARS_+4] = vz*nx - vx*nz;

	} else if (dir == _YDIR_) {
    un = vy;
    ny = 1.0;

		R[0*_MODEL_NVARS_+0] = 1.0;
		R[1*_MODEL_NVARS_+0] = vx - a*nx;
		R[2*_MODEL_NVARS_+0] = vy - a*ny;
		R[3*_MODEL_NVARS_+0] = vz - a*nz;
		R[4*_MODEL_NVARS_+0] = h0 - a*un;

		R[0*_MODEL_NVARS_+1] = 1.0;
		R[1*_MODEL_NVARS_+1] = vx;
		R[2*_MODEL_NVARS_+1] = vy;
		R[3*_MODEL_NVARS_+1] = vz;
		R[4*_MODEL_NVARS_+1] = ek;

		R[0*_MODEL_NVARS_+2] = 1.0;
		R[1*_MODEL_NVARS_+2] = vx + a*nx;
		R[2*_MODEL_NVARS_+2] = vy + a*ny;
		R[3*_MODEL_NVARS_+2] = vz + a*nz;
		R[4*_MODEL_NVARS_+2] = h0 + a*un;

		R[0*_MODEL_NVARS_+3] = 0;
		R[1*_MODEL_NVARS_+3] = ny;
		R[2*_MODEL_NVARS_+3] = -nx;
		R[3*_MODEL_NVARS_+3] = 0;
		R[4*_MODEL_NVARS_+3] =vx*ny-vy*nx;

		R[0*_MODEL_NVARS_+4] = 0;
		R[1*_MODEL_NVARS_+4] = 0;
		R[2*_MODEL_NVARS_+4] = nz;
		R[3*_MODEL_NVARS_+4] = -ny;
		R[4*_MODEL_NVARS_+4] = vy*nz-vz*ny;

  } else if (dir == _ZDIR_) { 
    un = vz; 
    nz = 1.0; 

		R[0*_MODEL_NVARS_+0] = 1.0;
		R[1*_MODEL_NVARS_+0] = vx-a*nx;
		R[2*_MODEL_NVARS_+0] = vy-a*ny;
		R[3*_MODEL_NVARS_+0] = vz-a*nz;
		R[4*_MODEL_NVARS_+0] = h0-a*un;

		R[0*_MODEL_NVARS_+1] = 1.0;
		R[1*_MODEL_NVARS_+1] = vx;
		R[2*_MODEL_NVARS_+1] = vy;
		R[3*_MODEL_NVARS_+1] = vz;
		R[4*_MODEL_NVARS_+1] = ek;

		R[0*_MODEL_NVARS_+2] = 1.0;
		R[1*_MODEL_NVARS_+2] = vx+a*nx;
		R[2*_MODEL_NVARS_+2] = vy+a*ny;
		R[3*_MODEL_NVARS_+2] = vz+a*nz;
		R[4*_MODEL_NVARS_+2] = h0+a*un;

		R[0*_MODEL_NVARS_+3] = 0;
		R[1*_MODEL_NVARS_+3] = -nz;
		R[2*_MODEL_NVARS_+3] = 0;
		R[3*_MODEL_NVARS_+3] = nx;
		R[4*_MODEL_NVARS_+3] = vz*nx-vx*nz;

		R[0*_MODEL_NVARS_+4] = 0;
		R[1*_MODEL_NVARS_+4] = 0;
		R[2*_MODEL_NVARS_+4] = nz;
		R[3*_MODEL_NVARS_+4] = -ny;
		R[4*_MODEL_NVARS_+4] = vy*nz-vz*ny;

  } else {
    fprintf(stderr,"Error in NavierStokes3DRightEigenvectors(): invalid dir!\n");
    return(1);
  }

  return(0);
}
