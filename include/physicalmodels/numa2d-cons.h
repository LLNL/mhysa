/*
 *
 * 2D Nonhydrostatic Unified Model of the Atmosphere (NUMA)
 * in the conservation form
 *
 * References:
 *
 * N. Ahmad and J. Lindeman, "Euler solutions using flux-based wave 
 * decomposition", Intl. J. Num. Method. Fluid., Vol. 54 (1), 
 * 47-72 (2007)
 *
*/

#include <basic.h>
#include <mathfunctions.h>

#define _NUMA2D_CONS_ "numa2d-cons"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 2
#define _MODEL_NVARS_ 4

/* grid directions */
#define _XDIR_ 0
#define _YDIR_ 1

#define _Numa2DConsComputePressure_(params,rho,theta,pressure) \
  { \
    double gamma    = params->gamma; \
    double Pref     = params->Pref; \
    double GasConst = params->R; \
    pressure  = Pref * raiseto((rho*GasConst*theta/Pref),gamma); \
  }

#define _Numa2DConsGetFlowVars_(u,rho,uvel,vvel,theta) \
  { \
    rho  = u[0]; \
    uvel = u[1] / rho; \
    vvel = u[2] / rho; \
    theta= u[3] / rho; \
  }

#define _Numa2DConsSetFlux_(f,dir,params,u) \
  { \
    double rho,uvel,vvel,theta,pressure; \
    _Numa2DConsGetFlowVars_(u,rho,uvel,vvel,theta); \
    _Numa2DConsComputePressure_(params,rho,theta,pressure); \
    if (dir == _XDIR_) { \
      f[0] = rho * uvel; \
      f[1] = rho*uvel*uvel + pressure; \
      f[2] = rho*uvel*vvel; \
      f[3] = rho * theta * uvel; \
    } else if (dir == _YDIR_) { \
      f[0] = rho * vvel; \
      f[1] = rho*uvel*vvel; \
      f[2] = rho*vvel*vvel + pressure; \
      f[3] = rho * theta * vvel; \
    } \
  }

#define _Numa2DConsSetLinearFlux_(f,dir,params,u) \
  { \
    double rho,uvel,vvel,theta,pressure; \
    _Numa2DConsGetFlowVars_(u,rho,uvel,vvel,theta); \
    _Numa2DConsComputePressure_(params,rho,theta,pressure); \
    if (dir == _XDIR_) { \
      f[0] = rho * uvel; \
      f[1] = pressure; \
      f[2] = 0.0; \
      f[3] = rho * theta * uvel; \
    } else if (dir == _YDIR_) { \
      f[0] = rho * vvel; \
      f[1] = 0.0; \
      f[2] = pressure; \
      f[3] = rho * theta * vvel; \
    } \
  }

#define _Numa2DConsSetSource_(s,param,u) \
  { \
    double rho,uvel,vvel,theta; \
    _Numa2DConsGetFlowVars_(u,rho,uvel,vvel,theta); \
    s[0] =  0.0; \
    s[1] =  0.0; \
    s[2] = -param->g*rho; \
    s[3] =  0.0; \
  }

/*
#define _Numa2DConsComputeLinearizedPressure_(params,T0,dT,P0,dP) \
  { \
    double gamma    = params->gamma; \
    dP  = (gamma*P0/T0) * dT; \
  }
*/

#define _Numa2DConsComputeSpeedofSound_(params,rho,theta,c) \
  { \
    double pressure, gamma = params->gamma; \
    _Numa2DConsComputePressure_(params,rho,theta,pressure); \
    c = sqrt(gamma*pressure/rho); \
  }
/*
#define _Numa2DConsRoeAverage_(uavg,u1,u2,params) \
  { \
    double gamma = params->gamma; \
    double GasConst = params->R; \
    double Pref     = params->Pref; \
    double rho1,uvel1,vvel1,theta1,c1,H1,t1; \
    _Numa2DConsGetFlowVars_(u1,rho1,uvel1,vvel1,theta1); \
    _Numa2DConsComputeSpeedofSound_(params,rho1,theta1,c1); \
    H1   = 0.5*(uvel1*uvel1+vvel1*vvel1) + c1*c1 / (gamma-1.0); \
    t1   = sqrt(rho1); \
    double rho2,uvel2,vvel2,theta2,c2,H2,t2; \
    _Numa2DConsGetFlowVars_(u2,rho2,uvel2,vvel2,theta2,); \
    _Numa2DConsComputeSpeedofSound_(params,rho2,theta2,c2); \
    H2   = 0.5*(uvel2*uvel2+vvel2*vvel2) + c2*c2 / (gamma-1.0); \
    t2   = sqrt(rho2); \
    double rho_avg,uvel_avg,vvel_avg,H_avg,c_sq_avg,pressure_avg,theta_avg; \
    rho_avg   = t1 * t2; \
    uvel_avg  = (t1*uvel1 + t2*uvel2) / (t1 + t2); \
    vvel_avg  = (t1*vvel1 + t2*vvel2) / (t1 + t2); \
    H_avg     = (t1*H1    + t2*H2   ) / (t1 + t2); \
    c_sq_avg = (gamma-1.0) * (H_avg - 0.5*(uvel_avg*uvel_avg+vvel_avg*vvel_avg)); \
    pressure_avg = c_sq_avg * rho_avg / gamma; \
    theta_avg = (Pref/(rho_avg*GasConst)) * raiseto((pressure_avg/Pref), (1.0/gamma)); \
    uavg[0] = rho_avg; \
    uavg[1] = rho_avg*uvel_avg; \
    uavg[2] = rho_avg*vvel_avg; \
    uavg[3] = rho_avg*theta_avg; \
  }
*/
typedef struct numa2dcons_parameters {
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
} Numa2DCons;

int Numa2DConsInitialize (void*,void*);
int Numa2DConsCleanup    (void*);

/* Available upwinding schemes */
#define _RUSANOV_UPWINDING_ "rusanov"
