/*
 *
 * 3D Nonhydrostatic Unified Model of the Atmosphere (NUMA)
 *
 * References:
 *
 * F.X. Giraldo, M. Restelli, and M. Laeuter, "Semi-Implicit 
 * Formulations of the Euler Equations: Applications to 
 * Nonhydrostatic Atmospheric Modeling", SIAM J. Sci. Comp., 
 * Vol. 32, 3394-3425 (2010)
 *
 * J.F. Kelly and F.X. Giraldo, "Continuous and Discontinuous 
 * Galerkin Methods for a Scalable 3D Nonhydrostatic Atmospheric 
 * Model:  limited-area mode", J. Comp. Phys., Vol. 231, 7988-8008 
 * (2012)
 *
 * N. Ahmad and J. Lindeman, "Euler solutions using flux-based wave 
 * decomposition", Intl. J. Num. Method. Fluid., Vol. 54 (1), 
 * 47-72 (2007)
 *
*/

#include <basic.h>
#include <mathfunctions.h>

#define _NUMA3D_ "numa3d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 3
#define _MODEL_NVARS_ 5

/* grid directions */
#define _XDIR_ 0
#define _YDIR_ 1
#define _ZDIR_ 2

#define _Numa3DGetFlowVars_(u,drho,uvel,vvel,wvel,dT,rho0) \
  { \
    drho = u[0]; \
    uvel = u[1]/(rho0+drho); \
    vvel = u[2]/(rho0+drho); \
    wvel = u[3]/(rho0+drho); \
    dT   = u[4]; \
  }

#define _Numa3DSetFlux_(f,dir,drho,uvel,vvel,wvel,dT,dP,rho0,T0) \
  { \
    if (dir == _XDIR_) { \
      f[0] = (rho0+drho) * uvel; \
      f[1] = (rho0+drho)*uvel*uvel + dP; \
      f[2] = (rho0+drho)*uvel*vvel; \
      f[3] = (rho0+drho)*uvel*wvel; \
      f[4] = uvel*(dT+T0); \
    } else if (dir == _YDIR_) { \
      f[0] = (rho0+drho) * vvel; \
      f[1] = (rho0+drho)*uvel*vvel; \
      f[2] = (rho0+drho)*vvel*vvel + dP; \
      f[3] = (rho0+drho)*wvel*vvel; \
      f[4] = vvel*(dT+T0); \
    } else if (dir == _ZDIR_) { \
      f[0] = (rho0+drho) * wvel; \
      f[1] = (rho0+drho)*uvel*wvel; \
      f[2] = (rho0+drho)*vvel*wvel; \
      f[3] = (rho0+drho)*wvel*wvel + dP; \
      f[4] = wvel*(dT+T0); \
    } \
  }

#define _Numa3DSetLinearFlux_(f,dir,drho,uvel,vvel,wvel,dT,dP,rho0,T0) \
  { \
    if (dir == _XDIR_) { \
      f[0] = (rho0+drho) * uvel; \
      f[1] = dP; \
      f[2] = 0.0; \
      f[3] = 0.0; \
      f[4] = (rho0+drho)*uvel*T0/rho0; \
    } else if (dir == _YDIR_) { \
      f[0] = (rho0+drho) * vvel; \
      f[1] = 0.0; \
      f[2] = dP; \
      f[3] = 0.0; \
      f[4] = (rho0+drho)*vvel*T0/rho0; \
    } else if (dir == _ZDIR_) { \
      f[0] = (rho0+drho) * wvel; \
      f[1] = 0.0; \
      f[2] = 0.0; \
      f[3] = dP; \
      f[4] = (rho0+drho)*wvel*T0/rho0; \
    } \
  }

#define _Numa3DSetSource_(s,param,uvel,vvel,drho,rho0) \
  { \
    s[0] =  0.0; \
    s[1] =  2.0*param->Omega*vvel*(rho0+drho); \
    s[2] = -2.0*param->Omega*uvel*(rho0+drho); \
    s[3] = -param->g*drho; \
    s[4] =  0.0; \
  }

#define _Numa3DComputePressure_(params,T0,dT,P0,dP) \
  { \
    double gamma    = params->gamma; \
    double Pref     = params->Pref; \
    double R        = params->R; \
    double P_total  = Pref * raiseto((R*(T0+dT)/Pref),gamma); \
    dP  = P_total - P0; \
  }

#define _Numa3DComputeLinearizedPressure_(params,T0,dT,P0,dP) \
  { \
    double gamma    = params->gamma; \
    dP  = (gamma*P0/T0) * dT; \
  }

#define _Numa3DComputeSpeedofSound_(gamma,R,T0,dT,rho0,drho,EP,c) \
  { \
    c = sqrt(gamma*R*(T0+dT)*EP/(rho0+drho)); \
  }

#define _Numa3DRoeAverage_(uavg,u1,u2,params,rho01,rho02,rho0,T01,T02,T0,EP1,EP2,EP) \
  { \
    double gamma = params->gamma; \
    double GasConst = params->R; \
    double drho1,uvel1,vvel1,wvel1,dT1,c1,rho1,H1,t1; \
    _Numa3DGetFlowVars_(u1,drho1,uvel1,vvel1,wvel1,dT1,rho01); \
    _Numa3DComputeSpeedofSound_(gamma,GasConst,T01,dT1,rho01,drho1,EP1,c1); \
    rho1 = rho01 + drho1; \
    H1   = 0.5*(uvel1*uvel1+vvel1*vvel1+wvel1*wvel1) + c1*c1 / (gamma-1.0); \
    t1   = sqrt(rho1); \
    double drho2,uvel2,vvel2,wvel2,dT2,c2,rho2,H2,t2; \
    _Numa3DGetFlowVars_(u2,drho2,uvel2,vvel2,wvel2,dT2,rho02); \
    _Numa3DComputeSpeedofSound_(gamma,GasConst,T02,dT2,rho02,drho2,EP2,c2); \
    rho2 = rho02 + drho2; \
    H2   = 0.5*(uvel2*uvel2+vvel2*vvel2+wvel2*wvel2) + c2*c2 / (gamma-1.0); \
    t2   = sqrt(rho2); \
    double rho_avg,uvel_avg,vvel_avg,wvel_avg,H_avg,c_sq_avg,T_avg; \
    rho_avg   = t1 * t2; \
    uvel_avg  = (t1*uvel1 + t2*uvel2) / (t1 + t2); \
    vvel_avg  = (t1*vvel1 + t2*vvel2) / (t1 + t2); \
    wvel_avg  = (t1*wvel1 + t2*wvel2) / (t1 + t2); \
    H_avg     = (t1*H1    + t2*H2   ) / (t1 + t2); \
    c_sq_avg = (gamma-1.0) * (H_avg - 0.5*(uvel_avg*uvel_avg+vvel_avg*vvel_avg+wvel_avg*wvel_avg)); \
    T_avg = ((c_sq_avg/(gamma*GasConst))/EP) * rho_avg; \
    uavg[0] = rho_avg-rho0; \
    uavg[1] = rho_avg*uvel_avg; \
    uavg[2] = rho_avg*vvel_avg; \
    uavg[3] = rho_avg*wvel_avg; \
    uavg[4] = T_avg-T0; \
  }

#define _Numa3DEigenvalues_(u,D,params,dir,rho0,T0,EP) \
  { \
    double gamma = params->gamma; \
    double GasConst = params->R; \
    double drho,uvel,vvel,wvel,dT,c,vn; \
    _Numa3DGetFlowVars_(u,drho,uvel,vvel,wvel,dT,rho0); \
    _Numa3DComputeSpeedofSound_(gamma,GasConst,T0,dT,rho0,drho,EP,c); \
    if      (dir == _XDIR_) vn = uvel; \
    else if (dir == _YDIR_) vn = vvel; \
    else if (dir == _ZDIR_) vn = wvel; \
    else                    vn = 0.0; \
    D[0*_MODEL_NVARS_+0] = vn-c; \
    D[0*_MODEL_NVARS_+1] = 0;    \
    D[0*_MODEL_NVARS_+2] = 0;    \
    D[0*_MODEL_NVARS_+3] = 0;    \
    D[0*_MODEL_NVARS_+4] = 0;    \
    D[1*_MODEL_NVARS_+0] = 0;    \
    D[1*_MODEL_NVARS_+1] = vn;   \
    D[1*_MODEL_NVARS_+2] = 0;    \
    D[1*_MODEL_NVARS_+3] = 0;    \
    D[1*_MODEL_NVARS_+4] = 0;    \
    D[2*_MODEL_NVARS_+0] = 0;    \
    D[2*_MODEL_NVARS_+1] = 0;    \
    D[2*_MODEL_NVARS_+2] = vn;   \
    D[2*_MODEL_NVARS_+3] = 0;    \
    D[2*_MODEL_NVARS_+4] = 0;    \
    D[3*_MODEL_NVARS_+0] = 0;    \
    D[3*_MODEL_NVARS_+1] = 0;    \
    D[3*_MODEL_NVARS_+2] = 0;    \
    D[3*_MODEL_NVARS_+3] = vn;   \
    D[3*_MODEL_NVARS_+4] = 0;    \
    D[4*_MODEL_NVARS_+0] = 0;    \
    D[4*_MODEL_NVARS_+1] = 0;    \
    D[4*_MODEL_NVARS_+2] = 0;    \
    D[4*_MODEL_NVARS_+3] = 0;    \
    D[4*_MODEL_NVARS_+4] = vn+c; \
  }

#define _Numa3DLeftEigenvectors_(u,L,params,dir,rho0,T0,EP) \
  { \
    double gamma    = params->gamma; \
    double GasConst = params->R; \
    double drho,uvel,vvel,wvel,dT,c; \
    _Numa3DGetFlowVars_(u,drho,uvel,vvel,wvel,dT,rho0); \
    _Numa3DComputeSpeedofSound_(gamma,GasConst,T0,dT,rho0,drho,EP,c); \
    double theta = (T0+dT) / (rho0+drho); \
    if (dir == _XDIR_) { \
    	L[0*_MODEL_NVARS_+0] = uvel/(2*c); \
  		L[0*_MODEL_NVARS_+1] = -1.0/(2*c); \
  		L[0*_MODEL_NVARS_+2] = 0; \
  		L[0*_MODEL_NVARS_+3] = 0; \
  		L[0*_MODEL_NVARS_+4] = 1.0/(2*theta); \
  		L[1*_MODEL_NVARS_+0] = 1.0; \
  		L[1*_MODEL_NVARS_+1] = 0; \
  		L[1*_MODEL_NVARS_+2] = 0; \
  		L[1*_MODEL_NVARS_+3] = 0; \
  		L[1*_MODEL_NVARS_+4] = -1.0/theta; \
  		L[2*_MODEL_NVARS_+0] = 0; \
  		L[2*_MODEL_NVARS_+1] = 0; \
  		L[2*_MODEL_NVARS_+2] = 1.0; \
  		L[2*_MODEL_NVARS_+3] = 0; \
  		L[2*_MODEL_NVARS_+4] = -vvel/theta; \
  		L[3*_MODEL_NVARS_+0] = 0; \
  		L[3*_MODEL_NVARS_+1] = 0; \
  		L[3*_MODEL_NVARS_+2] = 0; \
  		L[3*_MODEL_NVARS_+3] = 1.0; \
  		L[3*_MODEL_NVARS_+4] = -wvel/theta; \
  		L[4*_MODEL_NVARS_+0] = -uvel/(2*c); \
  		L[4*_MODEL_NVARS_+1] = 1.0/(2*c); \
  		L[4*_MODEL_NVARS_+2] = 0; \
  		L[4*_MODEL_NVARS_+3] = 0; \
  		L[4*_MODEL_NVARS_+4] = 1.0/(2*theta); \
    } else if (dir == _YDIR_) {  \
    	L[0*_MODEL_NVARS_+0] = vvel/(2*c); \
  		L[0*_MODEL_NVARS_+1] = 0; \
  		L[0*_MODEL_NVARS_+2] = -1.0/(2*c); \
  		L[0*_MODEL_NVARS_+3] = 0; \
  		L[0*_MODEL_NVARS_+4] = 1.0/(2*theta); \
  		L[1*_MODEL_NVARS_+0] = 1.0; \
  		L[1*_MODEL_NVARS_+1] = 0; \
  		L[1*_MODEL_NVARS_+2] = 0; \
  		L[1*_MODEL_NVARS_+3] = 0; \
  		L[1*_MODEL_NVARS_+4] = -1.0/theta; \
  		L[2*_MODEL_NVARS_+0] = 0; \
  		L[2*_MODEL_NVARS_+1] = 1.0; \
  		L[2*_MODEL_NVARS_+2] = 0; \
  		L[2*_MODEL_NVARS_+3] = 0; \
  		L[2*_MODEL_NVARS_+4] = -uvel/theta; \
  		L[3*_MODEL_NVARS_+0] = 0; \
  		L[3*_MODEL_NVARS_+1] = 0; \
  		L[3*_MODEL_NVARS_+2] = 0; \
  		L[3*_MODEL_NVARS_+3] = 1.0; \
  		L[3*_MODEL_NVARS_+4] = -wvel/theta; \
  		L[4*_MODEL_NVARS_+0] = -vvel/(2*c); \
  		L[4*_MODEL_NVARS_+1] = 0; \
  		L[4*_MODEL_NVARS_+2] = 1.0/(2*c); \
  		L[4*_MODEL_NVARS_+3] = 0; \
  		L[4*_MODEL_NVARS_+4] = 1.0/(2*theta); \
    } else if (dir == _ZDIR_) {  \
    	L[0*_MODEL_NVARS_+0] = wvel/(2*c); \
  		L[0*_MODEL_NVARS_+1] = 0; \
  		L[0*_MODEL_NVARS_+2] = 0; \
  		L[0*_MODEL_NVARS_+3] = -1.0/(2*c); \
  		L[0*_MODEL_NVARS_+4] = 1.0/(2*theta); \
  		L[1*_MODEL_NVARS_+0] = 1.0; \
  		L[1*_MODEL_NVARS_+1] = 0; \
  		L[1*_MODEL_NVARS_+2] = 0; \
  		L[1*_MODEL_NVARS_+3] = 0; \
  		L[1*_MODEL_NVARS_+4] = -1.0/theta; \
  		L[2*_MODEL_NVARS_+0] = 0; \
  		L[2*_MODEL_NVARS_+1] = 1.0; \
  		L[2*_MODEL_NVARS_+2] = 0; \
  		L[2*_MODEL_NVARS_+3] = 0; \
  		L[2*_MODEL_NVARS_+4] = -uvel/theta; \
  		L[3*_MODEL_NVARS_+0] = 0; \
  		L[3*_MODEL_NVARS_+1] = 0; \
  		L[3*_MODEL_NVARS_+2] = 1.0; \
  		L[3*_MODEL_NVARS_+3] = 0; \
  		L[3*_MODEL_NVARS_+4] = -vvel/theta; \
  		L[4*_MODEL_NVARS_+0] = -wvel/(2*c); \
  		L[4*_MODEL_NVARS_+1] = 0; \
  		L[4*_MODEL_NVARS_+2] = 0; \
  		L[4*_MODEL_NVARS_+3] = 1.0/(2*c); \
  		L[4*_MODEL_NVARS_+4] = 1.0/(2*theta); \
    } \
  }

#define _Numa3DRightEigenvectors_(u,R,params,dir,rho0,T0,EP) \
  { \
    double gamma    = params->gamma; \
    double GasConst = params->R; \
    double drho,uvel,vvel,wvel,dT,c; \
    _Numa3DGetFlowVars_(u,drho,uvel,vvel,wvel,dT,rho0); \
    _Numa3DComputeSpeedofSound_(gamma,GasConst,T0,dT,rho0,drho,EP,c); \
    double theta = (T0+dT) / (rho0+drho); \
	  if (dir == _XDIR_) { \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
		  R[1*_MODEL_NVARS_+0] = uvel-c; \
  		R[2*_MODEL_NVARS_+0] = vvel; \
	  	R[3*_MODEL_NVARS_+0] = wvel; \
		  R[4*_MODEL_NVARS_+0] = theta; \
  		R[0*_MODEL_NVARS_+1] = 1.0; \
	  	R[1*_MODEL_NVARS_+1] = uvel; \
		  R[2*_MODEL_NVARS_+1] = 0; \
  		R[3*_MODEL_NVARS_+1] = 0; \
	  	R[4*_MODEL_NVARS_+1] = 0; \
		  R[0*_MODEL_NVARS_+2] = 0; \
  		R[1*_MODEL_NVARS_+2] = 0; \
	  	R[2*_MODEL_NVARS_+2] = 1.0; \
		  R[3*_MODEL_NVARS_+2] = 0; \
  		R[4*_MODEL_NVARS_+2] = 0; \
	  	R[0*_MODEL_NVARS_+3] = 0; \
		  R[1*_MODEL_NVARS_+3] = 0; \
  		R[2*_MODEL_NVARS_+3] = 0; \
	  	R[3*_MODEL_NVARS_+3] = 1.0; \
		  R[4*_MODEL_NVARS_+3] = 0; \
  		R[0*_MODEL_NVARS_+4] = 1.0; \
	  	R[1*_MODEL_NVARS_+4] = uvel+c; \
		  R[2*_MODEL_NVARS_+4] = vvel; \
  		R[3*_MODEL_NVARS_+4] = wvel; \
	  	R[4*_MODEL_NVARS_+4] = theta; \
	  } else if (dir == _YDIR_) { \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
		  R[1*_MODEL_NVARS_+0] = uvel; \
  		R[2*_MODEL_NVARS_+0] = vvel-c; \
	  	R[3*_MODEL_NVARS_+0] = wvel; \
		  R[4*_MODEL_NVARS_+0] = theta; \
  		R[0*_MODEL_NVARS_+1] = 1.0; \
	  	R[1*_MODEL_NVARS_+1] = 0; \
		  R[2*_MODEL_NVARS_+1] = vvel; \
  		R[3*_MODEL_NVARS_+1] = 0; \
	  	R[4*_MODEL_NVARS_+1] = 0; \
		  R[0*_MODEL_NVARS_+2] = 0; \
  		R[1*_MODEL_NVARS_+2] = 1.0; \
	  	R[2*_MODEL_NVARS_+2] = 0; \
		  R[3*_MODEL_NVARS_+2] = 0; \
  		R[4*_MODEL_NVARS_+2] = 0; \
	  	R[0*_MODEL_NVARS_+3] = 0; \
		  R[1*_MODEL_NVARS_+3] = 0; \
  		R[2*_MODEL_NVARS_+3] = 0; \
	  	R[3*_MODEL_NVARS_+3] = 1.0; \
		  R[4*_MODEL_NVARS_+3] = 0; \
  		R[0*_MODEL_NVARS_+4] = 1.0; \
	  	R[1*_MODEL_NVARS_+4] = uvel; \
		  R[2*_MODEL_NVARS_+4] = vvel+c; \
  		R[3*_MODEL_NVARS_+4] = wvel; \
	  	R[4*_MODEL_NVARS_+4] = theta; \
    } else if (dir == _ZDIR_) {  \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
		  R[1*_MODEL_NVARS_+0] = uvel; \
  		R[2*_MODEL_NVARS_+0] = vvel; \
	  	R[3*_MODEL_NVARS_+0] = wvel-c; \
		  R[4*_MODEL_NVARS_+0] = theta; \
  		R[0*_MODEL_NVARS_+1] = 1.0; \
	  	R[1*_MODEL_NVARS_+1] = 0; \
		  R[2*_MODEL_NVARS_+1] = 0; \
  		R[3*_MODEL_NVARS_+1] = wvel; \
	  	R[4*_MODEL_NVARS_+1] = 0; \
		  R[0*_MODEL_NVARS_+2] = 0; \
  		R[1*_MODEL_NVARS_+2] = 1.0; \
	  	R[2*_MODEL_NVARS_+2] = 0; \
		  R[3*_MODEL_NVARS_+2] = 0; \
  		R[4*_MODEL_NVARS_+2] = 0; \
	  	R[0*_MODEL_NVARS_+3] = 0; \
		  R[1*_MODEL_NVARS_+3] = 0; \
  		R[2*_MODEL_NVARS_+3] = 1.0; \
	  	R[3*_MODEL_NVARS_+3] = 0; \
		  R[4*_MODEL_NVARS_+3] = 0; \
  		R[0*_MODEL_NVARS_+4] = 1.0; \
	  	R[1*_MODEL_NVARS_+4] = uvel; \
		  R[2*_MODEL_NVARS_+4] = vvel; \
  		R[3*_MODEL_NVARS_+4] = wvel+c; \
	  	R[4*_MODEL_NVARS_+4] = theta; \
    } \
  }

typedef struct numa3d_parameters {
  double  gamma;      /* Ratio of heat capacities       */
  double  R;          /* Universal gas constant         */
  double  Omega;      /* Angular speed of Earth         */
  double  g;          /* acceleration due to gravity    */
  int     init_atmos; /* choice of initial atmosphere   */

  /* pressure & temperature at zero altitude */
  double Pref, Tref;

  /* function to calculate hydrostatically balanced flow variables */
  void (*StandardAtmosphere)(void*,double,double*,double*,double*,double*);

  /* choice of upwinding scheme */
  char upwind[_MAX_STRING_SIZE_];
} Numa3D;

int Numa3DInitialize (void*,void*);
int Numa3DCleanup    (void*);

/* Available upwinding schemes */
#define _RUSANOV_UPWINDING_ "rusanov"
#define _RF_CHAR_UPWINDING_ "rf-char" /* does not work! */
