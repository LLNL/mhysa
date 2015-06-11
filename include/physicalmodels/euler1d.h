/*! @file euler1d.h
    @brief 1D Euler Equations (inviscid, compressible flows)
    @author Debojyoti Ghosh

  1D Euler Equations for Inviscid, Compressible Flows (includes gravitational force terms)\n

  \f{equation}{
    \frac {\partial} {\partial t} \left[\begin{array}{c} \rho \\ \rho u \\ e \end{array}\right]
  + \frac {\partial} {\partial x} \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ (e+p) u\end{array}\right]
  = \left[\begin{array}{c} 0 \\ -\rho g \\ -\rho u g \end{array}\right]
  \f}
  where
  \f{equation}{
    e = \frac {p} {\gamma-1} + \frac{1}{2} \rho u^2
  \f}

  For the treatment of gravitational source terms, refer to:\n
  + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
    for the Gas Dynamics Equations Under Gravitational Fields",
    Journal of Scientific Computing, 54, 2013, pp. 645-662.
    http://dx.doi.org/10.1007/s10915-012-9585-8
*/
/*
  d   [ rho   ]   d   [   rho*u    ]\n
  --  [ rho*u ] + --  [rho*u*u + p ] = 0\n
  dt  [   e   ]   dx  [ (e+p)*u    ]\n

  Equation of state:
           p         1
    e = -------  +   - rho * u^2
        gamma-1      2

*/

#include <basic.h>
#include <math.h>
#include <matops.h>

/*! \def _EULER_1D_ 
    1D Euler equations
*/
#define _EULER_1D_  "euler1d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 1
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 3

/* choices for upwinding schemes */
/*! Roe upwinding scheme */
#define _ROE_     "roe"
/*! Roe-Fixed upwinding scheme (characteristic-based Roe scheme with local Lax-Friedrich entropy fix) */
#define _RF_      "rf-char"
/*! Local Lax-Friedrich upwinding scheme */
#define _LLF_     "llf-char"
/*! Steger-Warming flux splitting */
#define _SWFS_    "steger-warming"

/* grid direction */
#undef _XDIR_
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0

/*! \def _Euler1DGetFlowVar_
  Compute the flow variables (density, pressure, velocity)
  from the conserved solution vector.
*/
#define _Euler1DGetFlowVar_(u,rho,v,e,P,p) \
  { \
    double gamma = p->gamma; \
    rho = u[0]; \
    v   = u[1] / rho; \
    e   = u[2]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
  }

/*! \def _Euler1DSetFlux_
    Set the flux vector given the flow variables (density,
    velocity, pressure).
*/
#define _Euler1DSetFlux_(f,rho,v,e,P) \
  { \
    f[0] = (rho) * (v); \
    f[1] = (rho) * (v) * (v) + (P); \
    f[2] = ((e) + (P)) * (v); \
  }

/*! \def _Euler1DSetStiffFlux_
    Set the stiff flux vector given the flow variables (density,
    velocity, pressure). The stiff flux vector is defines as the 
    part of the total flux that represents the two acoustic waves.
    This is used in a characteristic-based split-flux method for
    implicit-explicit time-integration. The stiff flux is integrated
    in time implicitly.
    \sa #_Euler1DSetStiffJac_, #_Euler1DSetLinearizedStiffFlux_
*/
#define _Euler1DSetStiffFlux_(f,rho,v,e,P,gamma) \
  { \
    double gamma_inv = 1.0 / (gamma); \
    f[0] = gamma_inv * (rho)*(v); \
    f[1] = gamma_inv * (rho)*(v)*(v) + (P); \
    f[2] = ((e)+(P)) * (v) - 0.5*gamma_inv*((gamma)-1)*(rho)*(v)*(v)*(v); \
  }

/*! \def _Euler1DSetLinearizedStiffFlux_
    Linearized version of #_Euler1DSetStiffFlux_, where the Jacobian
    of the stiff flux is evaluated at the beginning of each time
    step.
    \sa #_Euler1DSetStiffJac_
*/
#define _Euler1DSetLinearizedStiffFlux_(f,u,J) \
  { \
    _MatVecMultiply_((J),(u),(f),(_MODEL_NVARS_)); \
  }

/*! \def _Euler1DSetStiffJac_
    Compute the Jacobian of the stiff flux defined in
    #_Euler1DSetStiffFlux_. The Jacobian is stored in
    the row-major format.
*/
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

/*! \def _Euler1DRoeAverage_
    Compute the Roe average of two conserved solution
    vectors.
*/
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

/*! \def _Euler1DEigenvalues_
    Compute the eigenvalues, given the conserved solution vector. The
    eigenvalues are returned as a 3x3 matrix stored in row-major format.
    It is a diagonal matrix with the eigenvalues as diagonal elements. 
    Admittedly, it is a wasteful way of storing the eigenvalues.
*/
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

/*! \def _Euler1DLeftEigenvectors_
    Compute the matrix that has the left-eigenvectors as
    its rows. Stored in row-major format.
*/
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

/*! \def _Euler1DRightEigenvectors_
    Compute the matrix that has the right-eigenvectors as
    its columns. Stored in row-major format.
*/
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

/*! \def Euler1D
    \brief Structure containing variables and parameters specific to the 1D Euler equations.
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the 1D Euler equations.
*/

/*! \brief Structure containing variables and parameters specific to the 1D Euler equations.
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the 1D Euler equations.
*/
typedef struct euler1d_parameters {
  double  gamma;        /*!< Ratio of heat capacities (\f$\gamma\f$) */
  double  grav;         /*!< acceleration due to gravity */

  /*! type of gravitational field 
      (0 is isothermal with constant potential, 
       1 is isothermal with sinusoidal potential)\n
      See:
        + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
          for the Gas Dynamics Equations Under Gravitational Fields",
          Journal of Scientific Computing, 54, 2013, pp. 645-662.
          http://dx.doi.org/10.1007/s10915-012-9585-8
  */
  int     grav_type;    

  double  *grav_field;  /*!< Array to store the gravity potential field */
  double  *fast_jac;    /*!< Array to store the linearized fast-waves Jacobian over the whole domain */
  double  *solution;    /*!< Array to store the reference solution for linearization */
  char    upw_choice[_MAX_STRING_SIZE_]; /*!< Choice of upwinding scheme.\sa #_ROE_,#_LLF_,#_RF_,#_SWFS_ */

  /*! Function pointer to the function that computes the "upwinding" step in source term computation. To 
      understand the implementation of the gravitational source terms, see:
        + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
          for the Gas Dynamics Equations Under Gravitational Fields",
          Journal of Scientific Computing, 54, 2013, pp. 645-662.
          http://dx.doi.org/10.1007/s10915-012-9585-8
  */
  int     (*SourceUpwind)(double*,double*,double*,double*,int,void*,double);

} Euler1D;

/*! Function to initialize the 1D Euler module */
int    Euler1DInitialize (void*,void*);
/*! Function to clean up the 1D Euler module */
int    Euler1DCleanup    (void*);

