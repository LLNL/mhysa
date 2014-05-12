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

#define _Numa3DGetFlowVars_(u,drho,uvel,vvel,wvel,dT) \
  { \
    drho = u[0]; \
    uvel = u[1]; \
    vvel = u[2]; \
    wvel = u[3]; \
    dT   = u[4]; \
  }

#define _Numa3DSetFlux_(f,dir,drho,uvel,vvel,wvel,dT,dP,rho0) \
  { \
    if (dir == _XDIR_) { \
      f[0] = (drho+rho0) * uvel; \
      f[1] = uvel*uvel + dP/(drho+rho0); \
      f[2] = uvel*vvel; \
      f[3] = uvel*wvel; \
      f[4] = uvel*dT; \
    } else if (dir == _YDIR_) { \
      f[0] = (drho+rho0) * vvel; \
      f[1] = uvel*vvel; \
      f[2] = vvel*vvel + dP/(drho+rho0); \
      f[3] = wvel*vvel; \
      f[4] = vvel*dT; \
    } else if (dir == _ZDIR_) { \
      f[0] = (drho+rho0) * wvel; \
      f[1] = uvel*wvel; \
      f[2] = vvel*wvel; \
      f[3] = wvel*wvel + dP/(drho+rho0); \
      f[4] = wvel*dT; \
    } \
  }

#define _Numa3DSetSource_(s,param,uvel,vvel,drho,rho0) \
  { \
    s[0] =  0.0; \
    s[1] =  2.0*param->Omega*vvel; \
    s[2] = -2.0*param->Omega*uvel; \
    s[3] = -param->g*drho/(drho+rho0); \
    s[4] =  0.0; \
  }

#define _Numa3DComputePressure_(params,rho,T,P0,P,dP) \
  { \
    double gamma = params->gamma; \
    double Pref  = params->Pref; \
    double R     = params->R; \
    P   = Pref * raiseto((rho*R*T/Pref),gamma); \
    dP  = P - P0; \
  }

typedef struct numa3d_parameters {
  double  gamma;  /* Ratio of heat capacities                     */
  double  R;      /* Universal gas constant                       */
  double  Omega;  /* Angular speed of Earth                       */
  double  g;      /* acceleration due to gravity                  */

  /* pressure & temperature at zero altitude */
  double Pref, Tref;

  /* mean hydrostatic flow variables */
  double *rho0, *P0, *T0;
} Numa3D;

int Numa3DInitialize (void*,void*);
int Numa3DCleanup    (void*);
