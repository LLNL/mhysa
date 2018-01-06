/*! @file euler1d.h
    @brief 1D Multispecies Euler Equations (nonequilibrium reacting flows)
    @author Debojyoti Ghosh

  1D Multispecies Euler Equations for Nonequilibrium, Reacting Flows\n

  \f{equation}{
    \frac {\partial} {\partial t} \left[\begin{array}{c} \rho_0 \\ \vdots \\ \rho_{\left(n_s-1\right)} \\ \rho u \\ \rho E \\ \rho_0 E^v_0 \\ \vdots \\ \rho_{n_v} E^v_{n_v} \end{array}\right]
  + \frac {\partial} {\partial x} \left[\begin{array}{c} \rho_0 u \\ \vdots \\ \rho_{\left(n_s-1\right)} u \\ \rho u^2 + p \\ \left(\rho E + p\right) u \\ \rho_0 E^v_0 u \\ \vdots \\ \rho_{n_v} E^v_{n_v} u \end{array}\right]
  = \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ 0 \\ 0 \\ 0 \\ \vdots \\ 0 \end{array}\right]
  \f}
  where
  \f$n_s\f$ is the number of species, 
  \f$n_v\f$ is the number of vibrational energies to evolve, 
  \f$\rho_i\f$ are the densities of each species, 
  \f$E^v_i\f$ are the vibrational energies,
  \f$\rho = \sum_{i=0}^{\left(n_s-1\right)} \rho_i\f$ is the total density, and
  \f{equation}{
    E = \frac {1} {\gamma-1} \frac{p}{\rho} + \frac{1}{2} u^2
  \f}

  Reference for the governing equations:
  + Grossman, Cinnella, "Flux-Split Algorithms for Flows with Non-equilibrium Chemistry
    and Vibrational Relaxation", J. Comput. Phys., 88, pp. 131-168, 1990.
*/

#include <basic.h>
#include <math.h>
#include <matops.h>

/*! \def _EULER_1D_ 
    1D Multispecies Euler equations
*/
#define _EULER_1D_  "euler1d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 1

/* choices for upwind flux */
/*! Rusanov flux */
#define _RUSANOV_     "rusanov"

/* grid direction */
#undef _XDIR_
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0

/*! \def _Euler1DTotalDensity_
  Compute the total density from an array with species densities
*/
#define _Euler1DTotalDensity_(rho_t,rho_s,ns) \
  { \
    int mcounter; \
    rho_t = 0; \
    for (mcounter=0; mcounter<(ns); mcounter++) { \
      rho_t += rho_s[mcounter]; \
    } \
  }

/*! \def _Euler1DSpeedOfSound_
  Compute speed of sound from density and pressure
*/
#define _Euler1DSpeedOfSound_(c,gamma,P,rho) \
  { \
    c = sqrt(gamma*P/rho);\
  }

/*! \def _Euler1DComputePressure_
  Compute pressure from density, energy, and velocity
*/
#define _Euler1DComputePressure_(P,rho_s,rho_t,v,E,E_vib,gamma) \
  { \
    P = rho_t * (E-0.5*v*v) * (gamma-1.0); \
  }

/*! \def _Euler1DComputeTemperature_
  Compute temperature from density, velocity, pressure, and energy
*/
#define _Euler1DComputeTemperature_(T,rho_s,rho_t,v,P,E,E_vib,gamma) \
  { \
    T = P/rho_t; \
  }

/*! \def _Euler1DGetFlowVar_
  Compute the flow variables (density, pressure, velocity, temperature)
  from the conserved solution vector.
*/
#define _Euler1DGetFlowVar_(u,rho_s,rho_t,v,E,E_vib,P,T,ctxt) \
  { \
    int mcounter;\
    double gamma = ctxt->gamma; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      rho_s[mcounter] = u[mcounter]; \
    } \
    _Euler1DTotalDensity_(rho_t,rho_s,ns); \
    v = u[ns] / rho_t; \
    E = u[ns+1] / rho_t; \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      E_vib[mcounter] = u[ns+2+mcounter]/rho_s[mcounter]; \
    } \
    \
    _Euler1DComputePressure_(P,rho_s,rho_t,v,E,E_vib,gamma); \
    _Euler1DComputeTemperature_(T,rho_s,rho_t,v,P,E,E_vib,gamma); \
  }

/*! \def _Euler1DSetFlowVar_
  Set the conserved solution vector from the flow variables 
  (density, pressure, velocity)
*/
#define _Euler1DSetFlowVar_(u,rho_s,rho_t,v,E,E_vib,P,ctxt) \
  { \
    int mcounter;\
    double gamma = ctxt->gamma; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      u[mcounter] = rho_s[mcounter]; \
    } \
    u[ns] = rho_t*v; \
    u[ns+1] = rho_t*E; \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      u[ns+2+mcounter] = rho_s[mcounter]*E_vib[mcounter]; \
    } \
  }

/*! \def _Euler1DSetFlux_
    Set the flux vector given the flow variables (density,
    velocity, pressure).
*/
#define _Euler1DSetFlux_(f,rho_s,rho_t,v,E,E_vib,P,ctxt) \
  { \
    int mcounter; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      f[mcounter] = rho_s[mcounter] * (v); \
    } \
    f[ns] = (rho_t)*(v)*(v) + (P); \
    f[ns+1] = ((rho_t)*(E) + (P)) * (v); \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      f[ns+2+mcounter] = rho_s[mcounter]*E_vib[mcounter]*(v); \
    } \
  }

/*! \def _Euler1DSetSource_
  Set the source vector
*/
#define _Euler1DSetSource_(s,rhodot_s,Evdot,ctxt) \
  { \
    int mcounter; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      s[mcounter] = rhodot_s[mcounter]; \
    } \
    s[ns] = 0; \
    s[ns+1] = 0; \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      s[ns+2+mcounter] = Evdot[mcounter]; \
    } \
  }

/*! \def _Euler1DRoeAverage_
    Compute the Roe average of two conserved solution
    vectors.
*/
#define _Euler1DRoeAverage_(uavg,uL,uR,ctxt) \
  { \
    int mcounter; \
    double gamma = ctxt->gamma; \
    int    ns = ctxt->n_species; \
    int    nv = ctxt->n_vibeng; \
    \
    double rhoL[ns],rhoL_t,vL,EL,EL_vib[nv],PL,TL,HL,cLsq; \
    _Euler1DGetFlowVar_(uL,rhoL,rhoL_t,vL,EL,EL_vib,PL,TL,ctxt); \
    cLsq = gamma * PL/rhoL_t; \
    HL = 0.5*vL*vL + cLsq / (gamma-1.0); \
    double tL[ns], tL_t; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      tL[mcounter] = sqrt(rhoL[mcounter]); \
    } \
    tL_t = sqrt(rhoL_t); \
    \
    double rhoR[ns],rhoR_t,vR,ER,ER_vib[nv],PR,TR,HR,cRsq; \
    _Euler1DGetFlowVar_(uR,rhoR,rhoR_t,vR,ER,ER_vib,PR,TR,ctxt); \
    cRsq = gamma * PR/rhoR_t; \
    HR = 0.5*vR*vR + cRsq / (gamma-1.0); \
    double tR[ns], tR_t; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      tR[mcounter] = sqrt(rhoR[mcounter]); \
    } \
    tR_t = sqrt(rhoR_t); \
    \
    double rho_s[ns],rho_t,v,E,E_vib[nv],P,H,csq; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      rho_s[mcounter] = tL[mcounter] * tR[mcounter]; \
    } \
    rho_t = tL_t * tR_t; \
    v = (tL_t*vL + tR_t*vR) / (tL_t + tR_t); \
    H = (tL_t*HL + tR_t*HR) / (tL_t + tR_t); \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      E_vib[mcounter] = (tL_t*EL_vib[mcounter]+tR_t*ER_vib[mcounter])/(tL_t+tR_t); \
    } \
    csq = (gamma-1.0) * (H-0.5*v*v); \
    P   = csq * rho_t / gamma; \
    E   = 1/(gamma-1.0)*(P/rho_t) + 0.5*v*v; \
    \
    _Euler1DSetFlowVar_(uavg,rho_s,rho_t,v,E,E_vib,P,ctxt); \
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

  int     n_species;    /*!< Number of species */
  int     n_vibeng;     /*!< Number of vibrational energies */

  int     nvars;        /*!< Number of components of solution variables
                             (= #Euler1D::n_species+#Euler1D::n_vibeng+2), 
                             not an user input; it is set to #HyPar::nvars in 
                             Euler1DInitialize().*/

  double  gamma;        /*!< Ratio of heat capacities (\f$\gamma\f$) */

  char    upw_choice[_MAX_STRING_SIZE_]; /*!< Choice of upwinding scheme.\sa #_RUSANOV_ */
} Euler1D;

/*! Function to initialize the 1D Euler module */
int    Euler1DInitialize (void*,void*);
/*! Function to clean up the 1D Euler module */
int    Euler1DCleanup    (void*);

