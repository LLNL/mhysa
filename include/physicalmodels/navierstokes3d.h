/*! @file navierstokes3d.h
    @brief 3D Multispecies Navier Stokes equations (nonequilibrium reacting flows)
    @author Debojyoti Ghosh

  3D Multispecies Navier-Stokes Equations for Nonequilibrium, Reacting Flows\n

  \f{equation}{
    \frac {\partial} {\partial t} \left[\begin{array}{c} \rho_0 \\ \vdots \\ \rho_{n_s-1} \\ \rho u \\ \rho v \\ \rho w \\ \rho E \\ \rho_0 E^v_0 \\ \vdots \\ \rho_{n_v} E^v_{n_v} \end{array}\right]
  + \frac {\partial} {\partial x} \left[\begin{array}{c} \rho_0 u \\ \vdots \\ \rho_{n_s-1} u \\ \rho u^2 + p \\ \rho u v \\ \rho u w \\ (\rho E+p) u \\ \rho_0 E^v_0 u \\ \vdots \\ \rho_{n_v} E^v_{n_v} u \end{array}\right]
  + \frac {\partial} {\partial y} \left[\begin{array}{c} \rho_0 v \\ \vdots \\ \rho_{n_s-1} v \\ \rho u v \\ \rho v^2 + p \\ \rho v w \\ (\rho E+p) v \\ \rho_0 E^v_0 v \\ \vdots \\ \rho_{n_v} E^v_{n_v} v \end{array}\right]
  + \frac {\partial} {\partial z} \left[\begin{array}{c} \rho_0 w \\ \vdots \\ \rho_{n_s-1} w \\ \rho u w \\ \rho v w \\ \rho w^2 + p \\ (\rho E+p) w \\ \rho_0 E^v_0 w \\ \vdots \\ \rho_{n_v} E^v_{n_v} w \end{array}\right]
  = \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ \tau_{xx} \\ \tau_{yx} \\ \tau_{zx} \\ u \tau_{xx} + v \tau_{yx} + w \tau_{zx} - q_x \\ ?? \\ \vdots \\ ?? \end{array}\right]
  + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{zy} \\ u \tau_{xy} + v \tau_{yy} + w \tau_{zy} - q_y \\ ?? \\ \vdots \\ ??\end{array}\right]
  + \frac {\partial} {\partial z} \left[\begin{array}{c} 0 \\ \vdots \\ 0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\ u \tau_{xz} + v \tau_{yz} + w \tau_{zz} - q_z \\ ?? \\ \vdots \\ ??\end{array}\right]
  \f}
  where
  \f$n_s\f$ is the number of species, 
  \f$n_v\f$ is the number of vibrational energies to evolve, 
  \f$\rho_i\f$ are the densities of each species, 
  \f$E^v_i\f$ are the vibrational energies,
  \f$\rho = \sum_{i=0}^{\left(n_s-1\right)} \rho_i\f$ is the total density,
  \f${\bf \hat{i}},{\bf \hat{j}}, {\bf \hat{k}}\f$ are the unit vectors along the x, y, and z, the viscous terms are given by
  \f{align}{
    \tau_{ij} &= \frac{\mu}{Re_\infty} \left[ \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) - \frac{2}{3}\frac{\partial u_k}{\partial x_k} \delta_{ij} \right], \\
    q_i &= - \frac{\mu}{\left(\gamma-1\right)Re_\infty Pr} \frac{\partial T}{\partial x_i}
  \f}
  with \f$\mu\f$ being the viscosity coefficient (computed using Sutherland's law), and the equation of state is
  \f{equation}{
    E = \frac {1} {\gamma-1} \frac{p}{\rho} + \frac{1}{2} \left(u^2 + v^2 + w^2\right)
  \f}

  References for the governing equations (as well as non-dimensional form):-
  + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
    Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
    number is based on speed of sound, instead of the freestream velocity).

  More references need to be added here.

*/

#include <basic.h>

/*! 3D Navier-Stokes equations */
#define _NAVIER_STOKES_3D_  "navierstokes3d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 3

/* choices for upwinding schemes */
/*! Rusanov's upwinding scheme */
#define _RUSANOV_   "rusanov"

/* grid directions */
#undef _XDIR_
#undef _YDIR_
#undef _ZDIR_
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0
/*! dimension corresponding to the \a y spatial dimension */
#define _YDIR_ 1
/*! dimension corresponding to the \a z spatial dimension */
#define _ZDIR_ 2

/*! \def _NavierStokes3DTotalDensity_
  Compute the total density from an array with species densities
*/
#define _NavierStokes3DTotalDensity_(rho_t,rho_s,ns) \
  { \
    int mcounter; \
    rho_t = 0; \
    for (mcounter=0; mcounter<(ns); mcounter++) { \
      rho_t += rho_s[mcounter]; \
    } \
  }

/*! \def _NavierStokes3DSpeedOfSound_
  Compute speed of sound from density and pressure
*/
#define _NavierStokes3DSpeedOfSound_(c,gamma,P,rho) \
  { \
    c = sqrt((gamma)*(P)/(rho));\
  }

/*! \def _NavierStokes3DComputePressure_
  Compute pressure from density, energy, and velocity
*/
#define _NavierStokes3DComputePressure_(P,rho_s,rho_t,u,v,w,E,E_vib,gamma) \
  { \
    P = (rho_t) * ((E)-0.5*((u)*(u)+(v)*(v)+(w)*(w))) * ((gamma)-1.0); \
  }

/*! \def _NavierStokes3DComputeTemperature_
  Compute temperature from density, velocity, pressure, and energy
*/
#define _NavierStokes3DComputeTemperature_(T,rho_s,rho_t,u,v,w,P,E,E_vib,gamma) \
  { \
    T = (P)/(rho_t); \
  }

/*! \def _NavierStokes3DGetFlowVar_
 Get the flow variables from the conserved solution vector.
*/
#define _NavierStokes3DGetFlowVar_(u,rho_s,rho_t,vx,vy,vz,E,E_vib,P,T,ctxt) \
  { \
    int mcounter;\
    double gamma = ctxt->gamma; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      rho_s[mcounter] = u[mcounter]; \
    } \
    _NavierStokes3DTotalDensity_(rho_t,rho_s,ns); \
    vx = u[ns]  / rho_t; \
    vy = u[ns+1]/ rho_t; \
    vz = u[ns+2]/ rho_t; \
    E  = u[ns+3]/ rho_t; \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      E_vib[mcounter] = u[ns+4+mcounter]/rho_s[mcounter]; \
    } \
    \
    _NavierStokes3DComputePressure_(P,rho_s,rho_t,vx,vy,vz,E,E_vib,gamma); \
    _NavierStokes3DComputeTemperature_(T,rho_s,rho_t,vx,vy,vz,P,E,E_vib,gamma); \
  }

/*! \def _NavierStokes3DSetFlowVar_
  Set the conserved solution vector from the flow variables 
  (density, pressure, velocity)
*/
#define _NavierStokes3DSetFlowVar_(u,rho_s,rho_t,vx,vy,vz,E,E_vib,P,ctxt) \
  { \
    int mcounter;\
    double gamma = ctxt->gamma; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      (u)[mcounter] = (rho_s)[mcounter]; \
    } \
    (u)[ns]   = (rho_t)*(vx); \
    (u)[ns+1] = (rho_t)*(vy); \
    (u)[ns+2] = (rho_t)*(vz); \
    (u)[ns+3] = (rho_t)*(E); \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      (u)[ns+4+mcounter] = (rho_s)[mcounter]*(E_vib)[mcounter]; \
    } \
  }

/*! \def _NavierStokes3DSetFlux_
  Compute the flux vector, given the flow variables
*/
#define _NavierStokes3DSetFlux_(f,rho_s,rho_t,vx,vy,vz,E,E_vib,P,ctxt,dir) \
  { \
    int mcounter; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    if (dir == _XDIR_) { \
      for (mcounter = 0; mcounter < ns; mcounter++) { \
        f[mcounter] = rho_s[mcounter] * (vx); \
      } \
      f[ns]   = (rho_t)*(vx)*(vx) + (P); \
      f[ns+1] = (rho_t)*(vx)*(vy); \
      f[ns+2] = (rho_t)*(vx)*(vz); \
      f[ns+3] = ((rho_t)*(E) + (P)) * (vx); \
      for (mcounter = 0; mcounter < nv; mcounter++) { \
        f[ns+4+mcounter] = rho_s[mcounter]*E_vib[mcounter]*(vx); \
      } \
    } else if (dir == _YDIR_) { \
      for (mcounter = 0; mcounter < ns; mcounter++) { \
        f[mcounter] = rho_s[mcounter] * (vy); \
      } \
      f[ns]   = (rho_t)*(vy)*(vx); \
      f[ns+1] = (rho_t)*(vy)*(vy) + (P); \
      f[ns+2] = (rho_t)*(vy)*(vz); \
      f[ns+3] = ((rho_t)*(E) + (P)) * (vy); \
      for (mcounter = 0; mcounter < nv; mcounter++) { \
        f[ns+4+mcounter] = rho_s[mcounter]*E_vib[mcounter]*(vy); \
      } \
    } else if (dir == _ZDIR_) { \
      for (mcounter = 0; mcounter < ns; mcounter++) { \
        f[mcounter] = rho_s[mcounter] * (vz); \
      } \
      f[ns]   = (rho_t)*(vz)*(vx); \
      f[ns+1] = (rho_t)*(vz)*(vy); \
      f[ns+2] = (rho_t)*(vz)*(vz) + (P); \
      f[ns+3] = ((rho_t)*(E) + (P)) * (vz); \
      for (mcounter = 0; mcounter < nv; mcounter++) { \
        f[ns+4+mcounter] = rho_s[mcounter]*E_vib[mcounter]*(vz); \
      } \
    } \
  }

/*! \def _NavierStokes3DSetSource_
  Set the source vector
*/
#define _NavierStokes3DSetSource_(s,rhodot_s,Evdot,ctxt) \
  { \
    int mcounter; \
    int ns = ctxt->n_species; \
    int nv = ctxt->n_vibeng; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      s[mcounter] = rhodot_s[mcounter]; \
    } \
    s[ns]   = 0; \
    s[ns+1] = 0; \
    s[ns+2] = 0; \
    s[ns+3] = 0; \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      s[ns+4+mcounter] = Evdot[mcounter]; \
    } \
  }

/*! \def _NavierStokes3DRoeAverage_
  Compute the Roe-average of two solutions.
*/
#define _NavierStokes3DRoeAverage_(uavg,uL,uR,ctxt) \
  { \
    int mcounter; \
    double gamma = ctxt->gamma; \
    int    ns = ctxt->n_species; \
    int    nv = ctxt->n_vibeng; \
    \
    double rhoL[ns],rhoL_t,vxL,vyL,vzL,EL,EL_vib[nv],PL,TL,HL,cLsq; \
    _NavierStokes3DGetFlowVar_(uL,rhoL,rhoL_t,vxL,vyL,vzL,EL,EL_vib,PL,TL,ctxt); \
    cLsq = gamma * PL/rhoL_t; \
    HL = 0.5*(vxL*vxL+vyL*vyL+vzL*vzL) + cLsq / (gamma-1.0); \
    double tL[ns], tL_t; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      tL[mcounter] = sqrt(rhoL[mcounter]); \
    } \
    tL_t = sqrt(rhoL_t); \
    \
    double rhoR[ns],rhoR_t,vxR,vyR,vzR,ER,ER_vib[nv],PR,TR,HR,cRsq; \
    _NavierStokes3DGetFlowVar_(uR,rhoR,rhoR_t,vxR,vyR,vzR,ER,ER_vib,PR,TR,ctxt); \
    cRsq = gamma * PR/rhoR_t; \
    HR = 0.5*(vxR*vxR+vyR*vyR+vzR*vzR) + cRsq / (gamma-1.0); \
    double tR[ns], tR_t; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      tR[mcounter] = sqrt(rhoR[mcounter]); \
    } \
    tR_t = sqrt(rhoR_t); \
    \
    double rho_s[ns],rho_t,vx,vy,vz,E,E_vib[nv],P,H,csq; \
    for (mcounter = 0; mcounter < ns; mcounter++) { \
      rho_s[mcounter] = tL[mcounter] * tR[mcounter]; \
    } \
    rho_t = tL_t * tR_t; \
    vx = (tL_t*vxL + tR_t*vxR) / (tL_t + tR_t); \
    vy = (tL_t*vyL + tR_t*vyR) / (tL_t + tR_t); \
    vz = (tL_t*vzL + tR_t*vzR) / (tL_t + tR_t); \
    H = (tL_t*HL + tR_t*HR) / (tL_t + tR_t); \
    for (mcounter = 0; mcounter < nv; mcounter++) { \
      E_vib[mcounter] = (tL_t*EL_vib[mcounter]+tR_t*ER_vib[mcounter])/(tL_t+tR_t); \
    } \
    csq = (gamma-1.0) * (H-0.5*(vx*vx+vy*vy+vz*vz)); \
    P   = csq * rho_t / gamma; \
    E   = 1/(gamma-1.0)*(P/rho_t) + 0.5*(vx*vx+vy*vy+vz*vz); \
    \
    _NavierStokes3DSetFlowVar_(uavg,rho_s,rho_t,vx,vy,vz,E,E_vib,P,ctxt); \
  }

/*! \def NavierStokes3D
    \brief Structure containing variables and parameters specific to the 3D Navier Stokes equations.
 *  This structure contains the physical parameters, variables, and function pointers specific to 
 *  the 3D Navier-Stokes equations.
*/
/*! \brief Structure containing variables and parameters specific to the 3D Navier Stokes equations.
 *  This structure contains the physical parameters, variables, and function pointers specific to 
 *  the 3D Navier-Stokes equations.
*/
typedef struct navierstokes3d_parameters {

  int     n_species;    /*!< Number of species */
  int     n_vibeng;     /*!< Number of vibrational energies */

  int     nvars;        /*!< Number of components of solution variables
                             (= #NavierStokes3D::n_species+#NavierStokes3D::n_vibeng+4), 
                             not an user input; it is set to #HyPar::nvars in 
                             NavierStokes3DInitialize().*/

  double  gamma;  /*!< Ratio of heat capacities */
  double  Re;     /*!< Reynolds number */
  double  Pr;     /*!< Prandtl  number */
  double  Minf;   /*!< Freestream Mach number */
  double  C1,     /*!< Sutherlands law constant */
          C2;     /*!< Sutherlands law constant */

  char    upw_choice[_MAX_STRING_SIZE_]; /*!< choice of upwinding */

  char ib_write_surface_data[_MAX_STRING_SIZE_]; /*!< Flag to indicate whether to analyze and write surface data for
                                                      immersed body, if present. Applicable only if #HyPar::flag_ib is 1 */

} NavierStokes3D;

int    NavierStokes3DInitialize (void*,void*);
int    NavierStokes3DCleanup    (void*);

