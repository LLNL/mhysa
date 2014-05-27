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

#define _Numa3DComputeSpeedofSound_(gamma,R,T0,dT,rho0,drho,EP,c) \
  { \
    c = sqrt(gamma*R*(T0+dT)*EP/(rho0+drho)); \
  }

typedef struct numa3d_parameters {
  double  gamma;      /* Ratio of heat capacities       */
  double  R;          /* Universal gas constant         */
  double  Omega;      /* Angular speed of Earth         */
  double  g;          /* acceleration due to gravity    */
  int     init_atmos; /* choice of initial atmosphere   */

  /* pressure & temperature at zero altitude */
  double Pref, Tref;

  /* mean hydrostatic flow variables */
  double *rho0, *P0, *T0, *ExnerP;
} Numa3D;

int Numa3DInitialize (void*,void*);
int Numa3DCleanup    (void*);
