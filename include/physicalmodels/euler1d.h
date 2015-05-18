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


For the treatment of gravitational source terms, refer to:
+ Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
  for the Gas Dynamics Equations Under Gravitational Fields",
  Journal of Scientific Computing, 54, 2013, pp. 645-662
  http://dx.doi.org/10.1007/s10915-012-9585-8

*/

#include <basic.h>
#include <math.h>
#include <matops.h>

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

/* grid direction */
#undef _XDIR_
#define _XDIR_ 0

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

#define _Euler1DSetStiffFlux_(f,rho,v,e,P,gamma) \
  { \
    double gamma_inv = 1.0 / (gamma); \
    f[0] = gamma_inv * (rho)*(v); \
    f[1] = gamma_inv * (rho)*(v)*(v) + (P); \
    f[2] = ((e)+(P)) * (v) - 0.5*gamma_inv*((gamma)-1)*(rho)*(v)*(v)*(v); \
  }

#define _Euler1DSetLinearizedStiffFlux_(f,u,J) \
  { \
    _MatVecMultiply_((J),(u),(f),(_MODEL_NVARS_)); \
  }

/* row-major format */
#define _Euler1DSetStiffJac_(J,rho,v,e,P,gamma) \
  { \
    double gm1 = (gamma)-1; \
    double inv_gm = 1.0/(gamma); \
    double gm1_inv_gm = gm1 * inv_gm; \
    J[0] = 0.5*gm1_inv_gm*(rho)*(v)*(v)*(v)/(P) - (v); \
    J[1] = 1.0 - gm1_inv_gm*(rho)*(v)*(v)/(P); \
    J[2] = gm1_inv_gm*(rho)*(v)/(P); \
    J[3] = 0.5*gm1_inv_gm*(rho)*(v)*(v)*(v)*(v)/(P) + 0.5*((gamma)-5.0)*(v)*(v); \
    J[4] = (3.0-(gamma))*(v) - gm1_inv_gm*(rho)*(v)*(v)*(v)/(P); \
    J[5] = gm1_inv_gm*(rho)*(v)*(v)/(P) + gm1; \
    J[6] = 0.25*gm1_inv_gm*(rho)*(v)*(v)*(v)*(v)*(v)/(P) - (gamma)*(v)*(e)/(rho) + 0.5*(2.0*(gamma)-3.0)*(v)*(v)*(v); \
    J[7] = (gamma)*(e)/(rho) - 1.5*gm1*(v)*(v) - 0.5*gm1_inv_gm*(rho)*(v)*(v)*(v)*(v)/(P); \
    J[8] = 0.5*gm1_inv_gm*(rho)*(v)*(v)*(v)/(P) + (gamma)*(v); \
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
  double  grav;   /* acceleration due to gravity */
  int     grav_type; /* type of gravitational field */
  double  *grav_field; /* gravity potential field */
  double  *fast_jac; /* linearized fast-waves Jacobian */
  double  *solution; /* reference solution for linearization */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
  int     (*SourceUpwind)(double*,double*,double*,double*,int,void*,double);
} Euler1D;

int    Euler1DInitialize (void*,void*);
int    Euler1DCleanup    (void*);

