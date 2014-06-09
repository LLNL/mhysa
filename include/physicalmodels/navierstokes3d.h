/*

  3D Navier-Stokes Equations for Compressible Flows

  Choices for upwinding:
  "roe"       Roe upwinding
  "rf-char"   Roe-fixed
  "llf-char"  Local Lax-Friedrich

  Refer: Computational Fluid Mechanics and Heat Transfer
         by Tannehill, Anderson and Pletcher
         Chapter 5, Section 5.1.7 for the non-dimensional
         form of the NS equations.
*/

#include <basic.h>

#define _NAVIER_STOKES_3D_  "navierstokes3d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 3
#define _MODEL_NVARS_ 5

/* choices for upwinding schemes */
#define _ROE_   "roe"
#define _RF_    "rf-char"
#define _LLF_   "llf-char"

/* grid directions */
#define _XDIR_ 0
#define _YDIR_ 1
#define _ZDIR_ 2

#define _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p) \
  { \
    double gamma   = p->gamma, vsq; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    vz  = u[3] / rho; \
    e   = u[4]; \
    vsq  = vx*vx + vy*vy + vz*vz; \
    P   = (e - 0.5*rho*vsq) * (gamma-1.0); \
  }

#define _NavierStokes3DSetFlux_(f,rho,vx,vy,vz,e,P,dir) \
  { \
    if (dir == _XDIR_) { \
      f[0] = rho * vx; \
      f[1] = rho * vx * vx + P; \
      f[2] = rho * vx * vy; \
      f[3] = rho * vx * vz; \
      f[4] = (e + P) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = rho * vy; \
      f[1] = rho * vy * vx; \
      f[2] = rho * vy * vy + P; \
      f[3] = rho * vy * vz; \
      f[4] = (e + P) * vy; \
    } else if (dir == _ZDIR_) { \
      f[0] = rho * vz; \
      f[1] = rho * vz * vx; \
      f[2] = rho * vz * vy; \
      f[3] = rho * vz * vz + P; \
      f[4] = (e + P) * vz; \
    } \
  }

#define _NavierStokes3DRoeAverage_(uavg,uL,uR,p) \
  { \
    double  rho ,vx, vy, vz, e ,P ,H; \
    double  rhoL,vxL,vyL,vzL,eL,PL,HL,cLsq; \
    double  rhoR,vxR,vyR,vzR,eR,PR,HR,cRsq; \
    double  gamma = p->gamma; \
    _NavierStokes3DGetFlowVar_(uL,rhoL,vxL,vyL,vzL,eL,PL,p); \
    cLsq = gamma * PL/rhoL; \
    HL = 0.5*(vxL*vxL+vyL*vyL+vzL*vzL) + cLsq / (gamma-1.0); \
    _NavierStokes3DGetFlowVar_(uR,rhoR,vxR,vyR,vzR,eR,PR,p); \
    cRsq = gamma * PR/rhoR; \
    HR = 0.5*(vxR*vxR+vyR*vyR+vzR*vzR) + cRsq / (gamma-1.0); \
    double tL = sqrt(rhoL); \
    double tR = sqrt(rhoR); \
    rho = tL * tR; \
    vx  = (tL*vxL + tR*vxR) / (tL + tR); \
    vy  = (tL*vyL + tR*vyR) / (tL + tR); \
    vz  = (tL*vzL + tR*vzR) / (tL + tR); \
    H   = (tL*HL  + tR*HR ) / (tL + tR); \
    P = (H - 0.5* (vx*vx+vy*vy+vz*vz)) * (rho*(gamma-1.0))/gamma; \
    e   = P/(gamma-1.0) + 0.5*rho*(vx*vx+vy*vy+vz*vz); \
    uavg[0] = rho; \
    uavg[1] = rho*vx; \
    uavg[2] = rho*vy; \
    uavg[3] = rho*vz; \
    uavg[4] = e; \
  }

#define _NavierStokes3DEigenvalues_(u,D,p,dir) \
  { \
    double          gamma   = p->gamma; \
    double          rho,vx,vy,vz,e,P,c,vn; \
    _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p); \
    c    = sqrt(gamma*P/rho); \
    if      (dir == _XDIR_) vn = vx; \
    else if (dir == _YDIR_) vn = vy; \
    else if (dir == _ZDIR_) vn = vz; \
    else                    vn = 0.0; \
    D[0*_MODEL_NVARS_+0] = vn-c;   D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0;    D[0*_MODEL_NVARS_+4] = 0; \
    D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vn;   D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0;    D[1*_MODEL_NVARS_+4] = 0; \
    D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vn+c;   D[2*_MODEL_NVARS_+3] = 0;    D[2*_MODEL_NVARS_+4] = 0; \
    D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vn;   D[3*_MODEL_NVARS_+4] = 0; \
    D[4*_MODEL_NVARS_+0] = 0;      D[4*_MODEL_NVARS_+1] = 0;    D[4*_MODEL_NVARS_+2] = 0;      D[4*_MODEL_NVARS_+3] = 0;    D[4*_MODEL_NVARS_+4] = vn;\
  }

#define _NavierStokes3DLeftEigenvectors_(u,L,p,dir) \
  { \
    double  ga = p->gamma, ga_minus_one=ga-1.0; \
    double  rho,vx,vy,vz,e,P,a,un,ek; \
    double  nx = 0,ny = 0,nz = 0; \
    _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p); \
  	ek = 0.5 * (vx*vx + vy*vy + vz*vz); \
	  a = sqrt(ga * P / rho); \
    if (dir == _XDIR_) { \
      un = vx;  \
      nx = 1.0; \
    	L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek + a*un) / (2*a*a); \
  		L[0*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx - a*nx) / (2*a*a); \
  		L[0*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy - a*ny) / (2*a*a); \
  		L[0*_MODEL_NVARS_+3] = ((-ga_minus_one)*vz - a*nz) / (2*a*a); \
  		L[0*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[1*_MODEL_NVARS_+0] = (a*a - ga_minus_one*ek) / (a*a); \
  		L[1*_MODEL_NVARS_+1] = (ga_minus_one*vx) / (a*a); \
  		L[1*_MODEL_NVARS_+2] = (ga_minus_one*vy) / (a*a); \
  		L[1*_MODEL_NVARS_+3] = (ga_minus_one*vz) / (a*a); \
  		L[1*_MODEL_NVARS_+4] = (-ga_minus_one) / (a*a); \
  		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek - a*un) / (2*a*a); \
  		L[2*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx + a*nx) / (2*a*a); \
  		L[2*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy + a*ny) / (2*a*a); \
  		L[2*_MODEL_NVARS_+3] = ((-ga_minus_one)*vz + a*nz) / (2*a*a); \
  		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[3*_MODEL_NVARS_+0] = (vy - un*ny) / nx; \
  		L[3*_MODEL_NVARS_+1] = ny; \
  		L[3*_MODEL_NVARS_+2] = (ny*ny - 1.0) / nx; \
  		L[3*_MODEL_NVARS_+3] = ny * nz / nx; \
  		L[3*_MODEL_NVARS_+4] = 0.0; \
  		L[4*_MODEL_NVARS_+0] = (un*nz - vz) / nx; \
  		L[4*_MODEL_NVARS_+1] = - nz; \
  		L[4*_MODEL_NVARS_+2] = - ny * nz / nx; \
  		L[4*_MODEL_NVARS_+3] = (1 - nz*nz) / nx; \
  		L[4*_MODEL_NVARS_+4] = 0.0; \
    } else if (dir == _YDIR_) {  \
      un = vy;  \
      ny = 1.0;  \
  		L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*un) / (2*a*a); \
  		L[0*_MODEL_NVARS_+1] = ((1.0-ga)*vx - a*nx) / (2*a*a); \
  		L[0*_MODEL_NVARS_+2] = ((1.0-ga)*vy - a*ny) / (2*a*a); \
  		L[0*_MODEL_NVARS_+3] = ((1.0-ga)*vz - a*nz) / (2*a*a); \
  		L[0*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[1*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a); \
  		L[1*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a); \
  		L[1*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a); \
  		L[1*_MODEL_NVARS_+3] = ga_minus_one*vz / (a*a); \
  		L[1*_MODEL_NVARS_+4] = (1.0 - ga) / (a*a); \
  		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*un) / (2*a*a); \
  		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx + a*nx) / (2*a*a); \
  		L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy + a*ny) / (2*a*a); \
  		L[2*_MODEL_NVARS_+3] = ((1.0-ga)*vz + a*nz) / (2*a*a); \
  		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[3*_MODEL_NVARS_+0] = (un*nx-vx) / ny; \
  		L[3*_MODEL_NVARS_+1] = (1.0 - nx*nx) / ny; \
  		L[3*_MODEL_NVARS_+2] = - nx; \
  		L[3*_MODEL_NVARS_+3] = -nx*nz / ny; \
  		L[3*_MODEL_NVARS_+4] = 0; \
  		L[4*_MODEL_NVARS_+0] = (vz - un*nz) / ny; \
  		L[4*_MODEL_NVARS_+1] = nx*nz / ny; \
  		L[4*_MODEL_NVARS_+2] = nz; \
  		L[4*_MODEL_NVARS_+3] = (nz*nz - 1.0) / ny; \
  		L[4*_MODEL_NVARS_+4] = 0; \
    } else if (dir == _ZDIR_) {  \
      un = vz;  \
      nz = 1.0;  \
      L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*un) / (2*a*a); \
  		L[0*_MODEL_NVARS_+1] = ((1.0-ga)*vx-a*nx) / (2*a*a); \
  		L[0*_MODEL_NVARS_+2] = ((1.0-ga)*vy-a*ny) / (2*a*a); \
  		L[0*_MODEL_NVARS_+3] = ((1.0-ga)*vz-a*nz) / (2*a*a); \
  		L[0*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[1*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a); \
  		L[1*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a); \
  		L[1*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a); \
  		L[1*_MODEL_NVARS_+3] = ga_minus_one*vz / (a*a); \
  		L[1*_MODEL_NVARS_+4] = (1.0-ga) / (a*a); \
  		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*un) / (2*a*a); \
  		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx+a*nx) / (2*a*a); \
  		L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy+a*ny) / (2*a*a); \
  		L[2*_MODEL_NVARS_+3] = ((1.0-ga)*vz+a*nz) / (2*a*a); \
  		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[3*_MODEL_NVARS_+0] = (vx-un*nx) / nz; \
  		L[3*_MODEL_NVARS_+1] = (nx*nx-1.0) / nz; \
  		L[3*_MODEL_NVARS_+2] = nx*ny / nz; \
  		L[3*_MODEL_NVARS_+3] = nx; \
  		L[3*_MODEL_NVARS_+4] = 0; \
  		L[4*_MODEL_NVARS_+0] = (un*ny-vy) / nz; \
  		L[4*_MODEL_NVARS_+1] = -nx*ny/nz; \
  		L[4*_MODEL_NVARS_+2] = (1.0-ny*ny) / nz; \
  		L[4*_MODEL_NVARS_+3] = -ny; \
  		L[4*_MODEL_NVARS_+4] = 0; \
    } \
  }

#define _NavierStokes3DRightEigenvectors_(u,R,p,dir) \
  { \
    double  ga   = p->gamma, ga_minus_one = ga-1.0; \
    double  rho,vx,vy,vz,e,P,un,ek,a,h0; \
    double  nx = 0,ny = 0,nz = 0; \
    _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p); \
  	ek   = 0.5 * (vx*vx + vy*vy + vz*vz); \
	  a    = sqrt(ga * P / rho); \
    h0   = a*a / ga_minus_one + ek; \
	  if (dir == _XDIR_) { \
      un = vx; \
      nx = 1.0; \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
		  R[1*_MODEL_NVARS_+0] = vx - a*nx; \
  		R[2*_MODEL_NVARS_+0] = vy - a*ny; \
	  	R[3*_MODEL_NVARS_+0] = vz - a*nz; \
		  R[4*_MODEL_NVARS_+0] = h0 - a*un; \
  		R[0*_MODEL_NVARS_+1] = 1.0; \
	  	R[1*_MODEL_NVARS_+1] = vx; \
		  R[2*_MODEL_NVARS_+1] = vy; \
  		R[3*_MODEL_NVARS_+1] = vz; \
	  	R[4*_MODEL_NVARS_+1] = ek; \
		  R[0*_MODEL_NVARS_+2] = 1.0; \
  		R[1*_MODEL_NVARS_+2] = vx + a*nx; \
	  	R[2*_MODEL_NVARS_+2] = vy + a*ny; \
		  R[3*_MODEL_NVARS_+2] = vz + a*nz; \
  		R[4*_MODEL_NVARS_+2] = h0 + a*un; \
	  	R[0*_MODEL_NVARS_+3] = 0.0; \
		  R[1*_MODEL_NVARS_+3] = ny; \
  		R[2*_MODEL_NVARS_+3] = -nx; \
	  	R[3*_MODEL_NVARS_+3] = 0.0; \
		  R[4*_MODEL_NVARS_+3] = vx*ny - vy*nx; \
  		R[0*_MODEL_NVARS_+4] = 0.0; \
	  	R[1*_MODEL_NVARS_+4] = -nz; \
		  R[2*_MODEL_NVARS_+4] = 0.0; \
  		R[3*_MODEL_NVARS_+4] = nx; \
	  	R[4*_MODEL_NVARS_+4] = vz*nx - vx*nz; \
	  } else if (dir == _YDIR_) { \
      un = vy; \
      ny = 1.0; \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
  		R[1*_MODEL_NVARS_+0] = vx - a*nx; \
	  	R[2*_MODEL_NVARS_+0] = vy - a*ny; \
		  R[3*_MODEL_NVARS_+0] = vz - a*nz; \
  		R[4*_MODEL_NVARS_+0] = h0 - a*un; \
	  	R[0*_MODEL_NVARS_+1] = 1.0; \
  		R[1*_MODEL_NVARS_+1] = vx; \
	  	R[2*_MODEL_NVARS_+1] = vy; \
		  R[3*_MODEL_NVARS_+1] = vz; \
  		R[4*_MODEL_NVARS_+1] = ek; \
	  	R[0*_MODEL_NVARS_+2] = 1.0; \
		  R[1*_MODEL_NVARS_+2] = vx + a*nx; \
  		R[2*_MODEL_NVARS_+2] = vy + a*ny; \
	  	R[3*_MODEL_NVARS_+2] = vz + a*nz; \
		  R[4*_MODEL_NVARS_+2] = h0 + a*un; \
  		R[0*_MODEL_NVARS_+3] = 0; \
	  	R[1*_MODEL_NVARS_+3] = ny; \
		  R[2*_MODEL_NVARS_+3] = -nx; \
  		R[3*_MODEL_NVARS_+3] = 0; \
	  	R[4*_MODEL_NVARS_+3] =vx*ny-vy*nx; \
  		R[0*_MODEL_NVARS_+4] = 0; \
	  	R[1*_MODEL_NVARS_+4] = 0; \
		  R[2*_MODEL_NVARS_+4] = nz; \
  		R[3*_MODEL_NVARS_+4] = -ny; \
	  	R[4*_MODEL_NVARS_+4] = vy*nz-vz*ny; \
    } else if (dir == _ZDIR_) {  \
      un = vz;  \
      nz = 1.0;  \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
	  	R[1*_MODEL_NVARS_+0] = vx-a*nx; \
	  	R[2*_MODEL_NVARS_+0] = vy-a*ny; \
	  	R[3*_MODEL_NVARS_+0] = vz-a*nz; \
	  	R[4*_MODEL_NVARS_+0] = h0-a*un; \
	  	R[0*_MODEL_NVARS_+1] = 1.0; \
	  	R[1*_MODEL_NVARS_+1] = vx; \
	  	R[2*_MODEL_NVARS_+1] = vy; \
	  	R[3*_MODEL_NVARS_+1] = vz; \
	  	R[4*_MODEL_NVARS_+1] = ek; \
	  	R[0*_MODEL_NVARS_+2] = 1.0; \
	  	R[1*_MODEL_NVARS_+2] = vx+a*nx; \
	  	R[2*_MODEL_NVARS_+2] = vy+a*ny; \
	  	R[3*_MODEL_NVARS_+2] = vz+a*nz; \
	  	R[4*_MODEL_NVARS_+2] = h0+a*un; \
	  	R[0*_MODEL_NVARS_+3] = 0; \
	  	R[1*_MODEL_NVARS_+3] = -nz; \
	  	R[2*_MODEL_NVARS_+3] = 0; \
	  	R[3*_MODEL_NVARS_+3] = nx; \
	  	R[4*_MODEL_NVARS_+3] = vz*nx-vx*nz; \
	  	R[0*_MODEL_NVARS_+4] = 0; \
	  	R[1*_MODEL_NVARS_+4] = 0; \
	  	R[2*_MODEL_NVARS_+4] = nz; \
	  	R[3*_MODEL_NVARS_+4] = -ny; \
	  	R[4*_MODEL_NVARS_+4] = vy*nz-vz*ny; \
    } \
  }

typedef struct navierstokes3d_parameters {
  double  gamma;  /* Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
  double  Re;     /* Reynolds number */
  double  Pr;     /* Prandtl  number */
  double  Minf;   /* Freestream Mach number */
  double  C1,C2;  /* Sutherlands law constants */
} NavierStokes3D;

int    NavierStokes3DInitialize (void*,void*);
int    NavierStokes3DCleanup    (void*);

