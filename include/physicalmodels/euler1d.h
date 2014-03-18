/*

  1D Euler Equations for Inviscid, Compressible Flows

    
  d   [ rho   ]   d   [   rho*u    ]
  --  [ rho*u ] + --  [rho*u*u + p ] = 0
  dt  [   e   ]   dx  [ (e+p)*u    ]

  Equation of state:
           p         1
    e = -------  +   - rho * u^2
        gamma-1      2

  Choices for upwinding:
  "roe"       Roe upwinding
  "rf-char"   Roe-fixed
  "llf-char"  Local Lax-Friedrich

*/

#include <basic.h>
#include <math.h>

#define _EULER_1D_  "euler1d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 1
#define _MODEL_NVARS_ 3

/* choices for upwinding schemes */
#define _ROE_     "roe"
#define _RF_      "rf-char"
#define _LLF_     "llf-char"
#define _SWFS_    "steger-warming"

#define _Euler1DGetFlowVar_(u,rho,v,e,P,p) \
  { \
    double gamma = p->gamma; \
    rho = u[0]; \
    v   = u[1] / rho; \
    e   = u[2]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
  }

#define _Euler1DSetFlux_(f,rho,v,e,P) \
  { \
    f[0] = (rho) * (v); \
    f[1] = (rho) * (v) * (v) + (P); \
    f[2] = ((e) + (P)) * (v); \
  }

#define _Euler1DRoeAverage_(uavg,uL,uR,p) \
  { \
    double rho ,v ,e ,P ,H ,csq; \
    double rhoL,vL,eL,PL,HL,cLsq; \
    double rhoR,vR,eR,PR,HR,cRsq; \
    double gamma = p->gamma; \
    rhoL = uL[0]; \
    vL   = uL[1] / rhoL; \
    eL   = uL[2]; \
    PL   = (eL - 0.5*rhoL*vL*vL) * (gamma-1.0); \
    cLsq = gamma * PL/rhoL; \
    HL = 0.5*vL*vL + cLsq / (gamma-1.0); \
    rhoR = uR[0]; \
    vR   = uR[1] / rhoR; \
    eR   = uR[2]; \
    PR   = (eR - 0.5*rhoR*vR*vR) * (gamma-1.0); \
    cRsq = gamma * PR/rhoR; \
    HR = 0.5*vR*vR + cRsq / (gamma-1.0); \
    double tL = sqrt(rhoL); \
    double tR = sqrt(rhoR); \
    rho = tL * tR; \
    v   = (tL*vL + tR*vR) / (tL + tR); \
    H   = (tL*HL + tR*HR) / (tL + tR); \
    csq = (gamma-1.0) * (H-0.5*v*v); \
    P   = csq * rho / gamma; \
    e   = P/(gamma-1.0) + 0.5*rho*v*v; \
    \
    uavg[0] = rho; \
    uavg[1] = rho*v; \
    uavg[2] = e; \
  }

#define _Euler1DEigenvalues_(u,D,p,dir) \
  { \
    double gamma   = p->gamma; \
    double rho,v,e,P,c; \
    rho = u[0]; \
    v   = u[1] / rho; \
    e   = u[2]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
    c = sqrt(gamma*P/rho); \
    D[0*_MODEL_NVARS_+0] = v;      D[0*_MODEL_NVARS_+1] = 0;      D[0*_MODEL_NVARS_+2] = 0; \
    D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = (v-c);  D[1*_MODEL_NVARS_+2] = 0; \
    D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;      D[2*_MODEL_NVARS_+2] = (v+c); \
  }

#define _Euler1DLeftEigenvectors_(u,L,p,dir) \
  { \
    double gamma = p->gamma; \
    double rho,v,e,P,c; \
    rho = u[0]; \
    v   = u[1] / rho; \
    e   = u[2]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
    c    = sqrt(gamma*P/rho); \
    L[1*_MODEL_NVARS_+0] = ((gamma - 1)/(rho*c)) * (-(v*v)/2 - c*v/(gamma-1)); \
    L[1*_MODEL_NVARS_+1] = ((gamma - 1)/(rho*c)) * (v + c/(gamma-1)); \
    L[1*_MODEL_NVARS_+2] = ((gamma - 1)/(rho*c)) * (-1); \
    L[0*_MODEL_NVARS_+0] = ((gamma - 1)/(rho*c)) * (rho*(-(v*v)/2+c*c/(gamma-1))/c); \
    L[0*_MODEL_NVARS_+1] = ((gamma - 1)/(rho*c)) * (rho*v/c); \
    L[0*_MODEL_NVARS_+2] = ((gamma - 1)/(rho*c)) * (-rho/c); \
    L[2*_MODEL_NVARS_+0] = ((gamma - 1)/(rho*c)) * ((v*v)/2 - c*v/(gamma-1)); \
    L[2*_MODEL_NVARS_+1] = ((gamma - 1)/(rho*c)) * (-v + c/(gamma-1)); \
    L[2*_MODEL_NVARS_+2] = ((gamma - 1)/(rho*c)) * (1); \
  }

#define _Euler1DRightEigenvectors_(u,R,p,dir) \
  { \
    double gamma   = p->gamma; \
    double rho,v,e,P,c; \
    rho = u[0]; \
    v   = u[1] / rho; \
    e   = u[2]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
    c    = sqrt(gamma*P/rho); \
    R[0*_MODEL_NVARS_+1] = - rho/(2*c);  R[1*_MODEL_NVARS_+1] = -rho*(v-c)/(2*c); R[2*_MODEL_NVARS_+1] = -rho*((v*v)/2+(c*c)/(gamma-1)-c*v)/(2*c);  \
    R[0*_MODEL_NVARS_+0] = 1;            R[1*_MODEL_NVARS_+0] = v;                R[2*_MODEL_NVARS_+0] = v*v / 2;                                   \
    R[0*_MODEL_NVARS_+2] = rho/(2*c);    R[1*_MODEL_NVARS_+2] = rho*(v+c)/(2*c);  R[2*_MODEL_NVARS_+2] = rho*((v*v)/2+(c*c)/(gamma-1)+c*v)/(2*c);   \
  }

typedef struct euler1d_parameters {
  double  gamma;  /* Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
} Euler1D;

int    Euler1DInitialize (void*,void*);
int    Euler1DCleanup    (void*);

