/*

  2D Navier-Stokes Equations for Inviscid, Compressible Flows

  This is a copy of the Euler2D physical system, along with 
  viscous terms.

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

#define _NAVIER_STOKES_2D_  "navierstokes2d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 2
#define _MODEL_NVARS_ 4

/* choices for upwinding schemes */
#define _ROE_   "roe"
#define _RF_    "rf-char"
#define _LLF_   "llf-char"
#define _SWFS_  "steger-warming"

/* directions */
#define _XDIR_ 0
#define _YDIR_ 1


#define _NavierStokes2DGetFlowVar_(u,rho,vx,vy,e,P,p) \
  { \
    double  gamma = p->gamma, vsq; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (gamma-1.0); \
  }

#define _NavierStokes2DSetFlux_(f,rho,vx,vy,e,P,p,dir) \
  { \
    if (dir == _XDIR_) { \
      f[0] = rho * vx; \
      f[1] = rho * vx * vx + P; \
      f[2] = rho * vx * vy; \
      f[3] = (e + P) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = rho * vy; \
      f[1] = rho * vy * vx; \
      f[2] = rho * vy * vy + P; \
      f[3] = (e + P) * vy; \
    } \
  }

#define _NavierStokes2DRoeAverage_(uavg,uL,uR,p) \
  { \
    double  rho ,vx, vy, e ,P ,H ,csq, vsq; \
    double  rhoL,vxL,vyL,eL,PL,HL,cLsq,vsqL; \
    double  rhoR,vxR,vyR,eR,PR,HR,cRsq,vsqR; \
    double  gamma = p->gamma; \
    rhoL = uL[0]; \
    vxL  = uL[1] / rhoL; \
    vyL  = uL[2] / rhoL; \
    eL   = uL[3]; \
    vsqL = (vxL*vxL) + (vyL*vyL); \
    PL   = (eL - 0.5*rhoL*vsqL) * (gamma-1.0); \
    cLsq = gamma * PL/rhoL; \
    HL = 0.5*(vxL*vxL+vyL*vyL) + cLsq / (gamma-1.0); \
    rhoR = uR[0]; \
    vxR  = uR[1] / rhoR; \
    vyR  = uR[2] / rhoR; \
    eR   = uR[3]; \
    vsqR = (vxR*vxR) + (vyR*vyR); \
    PR   = (eR - 0.5*rhoR*vsqR) * (gamma-1.0); \
    cRsq = gamma * PR/rhoR; \
    HR = 0.5*(vxR*vxR+vyR*vyR) + cRsq / (gamma-1.0); \
    double tL = sqrt(rhoL); \
    double tR = sqrt(rhoR); \
    rho = tL * tR; \
    vx  = (tL*vxL + tR*vxR) / (tL + tR); \
    vy  = (tL*vyL + tR*vyR) / (tL + tR); \
    H   = (tL*HL + tR*HR) / (tL + tR); \
    vsq = vx*vx + vy*vy; \
    csq = (gamma-1.0) * (H-0.5*vsq); \
    P   = csq * rho / gamma; \
    e   = P/(gamma-1.0) + 0.5*rho*vsq; \
    uavg[0] = rho; \
    uavg[1] = rho*vx; \
    uavg[2] = rho*vy; \
    uavg[3] = e; \
  }

#define _NavierStokes2DEigenvalues_(u,D,p,dir) \
  { \
    double  gamma = p->gamma; \
    double  rho,vx,vy,e,P,c,vn,vsq; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (gamma-1.0); \
    c    = sqrt(gamma*P/rho); \
    if      (dir == _XDIR_) vn = vx; \
    else if (dir == _YDIR_) vn = vy; \
    else               vn = 0.0; \
    if (dir == _XDIR_) {\
      D[0*_MODEL_NVARS_+0] = vn-c;   D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vn+c; D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vn;     D[2*_MODEL_NVARS_+3] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vn; \
    } else if (dir == _YDIR_) { \
      D[0*_MODEL_NVARS_+0] = vn-c;   D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vn;   D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vn+c;   D[2*_MODEL_NVARS_+3] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vn; \
    }\
  }

#define _NavierStokes2DLeftEigenvectors_(u,L,p,dir) \
  { \
    double  ga = param->gamma, ga_minus_one=ga-1.0; \
    double  rho,vx,vy,e,P,a,un,ek,vsq; \
    double  nx = 0,ny = 0; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (ga-1.0); \
  	ek = 0.5 * (vx*vx + vy*vy); \
	  a = sqrt(ga * P / rho); \
    if (dir == _XDIR_) { \
      un = vx; \
      nx = 1.0; \
  		L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek + a*un) / (2*a*a); \
	  	L[0*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx - a*nx) / (2*a*a); \
		  L[0*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy - a*ny) / (2*a*a); \
  		L[0*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
	  	L[3*_MODEL_NVARS_+0] = (a*a - ga_minus_one*ek) / (a*a); \
  		L[3*_MODEL_NVARS_+1] = (ga_minus_one*vx) / (a*a); \
	  	L[3*_MODEL_NVARS_+2] = (ga_minus_one*vy) / (a*a); \
		  L[3*_MODEL_NVARS_+3] = (-ga_minus_one) / (a*a); \
		  L[1*_MODEL_NVARS_+0] = (ga_minus_one*ek - a*un) / (2*a*a); \
  		L[1*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx + a*nx) / (2*a*a); \
	  	L[1*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy + a*ny) / (2*a*a); \
		  L[1*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
  		L[2*_MODEL_NVARS_+0] = (vy - un*ny) / nx; \
	  	L[2*_MODEL_NVARS_+1] = ny; \
		  L[2*_MODEL_NVARS_+2] = (ny*ny - 1.0) / nx; \
  		L[2*_MODEL_NVARS_+3] = 0.0; \
    } else if (dir == _YDIR_) {  \
      un = vy;  \
      ny = 1.0; \
	  	L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*un) / (2*a*a); \
		  L[0*_MODEL_NVARS_+1] = ((1.0-ga)*vx - a*nx) / (2*a*a); \
  		L[0*_MODEL_NVARS_+2] = ((1.0-ga)*vy - a*ny) / (2*a*a); \
	  	L[0*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
		  L[3*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a); \
  		L[3*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a); \
	  	L[3*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a); \
		  L[3*_MODEL_NVARS_+3] = (1.0 - ga) / (a*a); \
		  L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*un) / (2*a*a); \
  		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx + a*nx) / (2*a*a); \
	  	L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy + a*ny) / (2*a*a); \
		  L[2*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
		  L[1*_MODEL_NVARS_+0] = (un*nx-vx) / ny; \
  		L[1*_MODEL_NVARS_+1] = (1.0 - nx*nx) / ny; \
	  	L[1*_MODEL_NVARS_+2] = - nx; \
		  L[1*_MODEL_NVARS_+3] = 0; \
    } \
  }

#define _NavierStokes2DRightEigenvectors_(u,R,p,dir) \
  { \
    double  ga   = param->gamma, ga_minus_one = ga-1.0; \
    double  rho,vx,vy,e,P,un,ek,a,h0,vsq; \
    double  nx = 0,ny = 0; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (ga-1.0); \
	  ek   = 0.5 * (vx*vx + vy*vy); \
  	a    = sqrt(ga * P / rho); \
    h0   = a*a / ga_minus_one + ek; \
	  if (dir == _XDIR_) { \
    	un = vx; \
      nx = 1.0; \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
	  	R[1*_MODEL_NVARS_+0] = vx - a*nx; \
	  	R[2*_MODEL_NVARS_+0] = vy - a*ny; \
	  	R[3*_MODEL_NVARS_+0] = h0 - a*un; \
  		R[0*_MODEL_NVARS_+3] = 1.0; \
  		R[1*_MODEL_NVARS_+3] = vx; \
  		R[2*_MODEL_NVARS_+3] = vy; \
  		R[3*_MODEL_NVARS_+3] = ek; \
  		R[0*_MODEL_NVARS_+1] = 1.0; \
  		R[1*_MODEL_NVARS_+1] = vx + a*nx; \
  		R[2*_MODEL_NVARS_+1] = vy + a*ny; \
  		R[3*_MODEL_NVARS_+1] = h0 + a*un; \
  		R[0*_MODEL_NVARS_+2] = 0.0; \
  		R[1*_MODEL_NVARS_+2] = ny; \
  		R[2*_MODEL_NVARS_+2] = -nx; \
  		R[3*_MODEL_NVARS_+2] = vx*ny - vy*nx; \
  	} else if (dir == _YDIR_) { \
      un = vy; \
      ny = 1.0; \
  		R[0*_MODEL_NVARS_+0] = 1.0; \
  		R[1*_MODEL_NVARS_+0] = vx - a*nx; \
  		R[2*_MODEL_NVARS_+0] = vy - a*ny; \
  		R[3*_MODEL_NVARS_+0] = h0 - a*un; \
  		R[0*_MODEL_NVARS_+3] = 1.0; \
  		R[1*_MODEL_NVARS_+3] = vx; \
  		R[2*_MODEL_NVARS_+3] = vy; \
  		R[3*_MODEL_NVARS_+3] = ek; \
  		R[0*_MODEL_NVARS_+2] = 1.0; \
  		R[1*_MODEL_NVARS_+2] = vx + a*nx; \
  		R[2*_MODEL_NVARS_+2] = vy + a*ny; \
  		R[3*_MODEL_NVARS_+2] = h0 + a*un; \
  		R[0*_MODEL_NVARS_+1] = 0; \
  		R[1*_MODEL_NVARS_+1] = ny; \
  		R[2*_MODEL_NVARS_+1] = -nx; \
  		R[3*_MODEL_NVARS_+1] = vx*ny-vy*nx; \
    } \
  }

typedef struct navierstokes2d_parameters {
  double  gamma;  /* Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
  double  Re;     /* Reynolds number */
  double  Pr;     /* Prandtl  number */
  double  Minf;   /* Freestream Mach number */
  double  C1,C2;  /* Sutherlands law constants */
} NavierStokes2D;

int    NavierStokes2DInitialize (void*,void*);
int    NavierStokes2DCleanup    (void*);

