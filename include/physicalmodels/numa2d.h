/*
 *
 * 2D Nonhydrostatic Unified Model of the Atmosphere (NUMA)
 *
 * References:
 *
 * F.X. Giraldo, M. Restelli, and M. Laeuter, "Semi-Implicit 
 * Formulations of the Euler Equations: Applications to 
 * Nonhydrostatic Atmospheric Modeling", SIAM J. Sci. Comp., 
 * Vol. 32, 3394-3425 (2010)
 *
 * N. Ahmad and J. Lindeman, "Euler solutions using flux-based wave 
 * decomposition", Intl. J. Num. Method. Fluid., Vol. 54 (1), 
 * 47-72 (2007)
 *
*/

#include <basic.h>
#include <mathfunctions.h>

#define _NUMA2D_ "numa2d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 2
#define _MODEL_NVARS_ 4

/* grid directions */
#define _XDIR_ 0
#define _YDIR_ 1

#define _Numa2DGetFlowVars_(u,drho,uvel,vvel,dT,rho0) \
  { \
    drho = u[0]; \
    uvel = u[1]/(rho0+drho); \
    vvel = u[2]/(rho0+drho); \
    dT   = u[3]; \
  }

#define _Numa2DSetFlux_(f,dir,drho,uvel,vvel,dT,dP,rho0,T0) \
  { \
    if (dir == _XDIR_) { \
      f[0] = (rho0+drho) * uvel; \
      f[1] = (rho0+drho)*uvel*uvel + dP; \
      f[2] = (rho0+drho)*uvel*vvel; \
      f[3] = uvel*(dT+T0); \
    } else if (dir == _YDIR_) { \
      f[0] = (rho0+drho) * vvel; \
      f[1] = (rho0+drho)*uvel*vvel; \
      f[2] = (rho0+drho)*vvel*vvel + dP; \
      f[3] = vvel*(dT+T0); \
    } \
  }

#define _Numa2DSetLinearFlux_(f,dir,drho,uvel,vvel,dT,dP,rho0,T0) \
  { \
    if (dir == _XDIR_) { \
      f[0] = (rho0+drho) * uvel; \
      f[1] = dP; \
      f[2] = 0.0; \
      f[3] = (rho0+drho)*uvel*T0/rho0; \
    } else if (dir == _YDIR_) { \
      f[0] = (rho0+drho) * vvel; \
      f[1] = 0.0; \
      f[2] = dP; \
      f[3] = (rho0+drho)*vvel*T0/rho0; \
    } \
  }

#define _Numa2DSetSource_(s,param,drho) \
  { \
    s[0] =  0.0; \
    s[1] =  0.0; \
    s[2] = -param->g*drho; \
    s[3] =  0.0; \
  }

#define _Numa2DComputePressure_(params,T0,dT,P0,dP) \
  { \
    double gamma    = params->gamma; \
    double Pref     = params->Pref; \
    double R        = params->R; \
    double P_total  = Pref * raiseto((R*(T0+dT)/Pref),gamma); \
    dP  = P_total - P0; \
  }

#define _Numa2DComputeLinearizedPressure_(params,T0,dT,P0,dP) \
  { \
    double gamma    = params->gamma; \
    dP  = (gamma*P0/T0) * dT; \
  }

#define _Numa2DComputeSpeedofSound_(gamma,R,T0,dT,rho0,drho,EP,c) \
  { \
    c = sqrt(gamma*R*(T0+dT)*EP/(rho0+drho)); \
  }

#define _Numa2DComputeLinearizedSpeedofSound_(gamma,R,T0,rho0,EP,c) \
  { \
    c = sqrt(gamma*R*T0*EP/rho0); \
  }

typedef struct numa2d_parameters {
  double  gamma;      /* Ratio of heat capacities       */
  double  R;          /* Universal gas constant         */
  double  g;          /* acceleration due to gravity    */
  int     init_atmos; /* choice of initial atmosphere   */

  /* pressure & temperature at zero altitude */
  double Pref, Tref;

  /* function to calculate hydrostatically balanced flow variables */
  void (*StandardAtmosphere)(void*,double,double*,double*,double*,double*);

  /* choice of upwinding scheme */
  char upwind[_MAX_STRING_SIZE_];
} Numa2D;

int Numa2DInitialize (void*,void*);
int Numa2DCleanup    (void*);

/* Available upwinding schemes */
#define _RUSANOV_UPWINDING_ "rusanov"
