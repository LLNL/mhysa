/*! @file navierstokes3d.h
    @brief 3D Navier Stokes equations (compressible flows)
    @author Debojyoti Ghosh

  3D Navier-Stokes equations for viscous and inviscid compressible flows (with gravitational terms)\n

  \f{equation}{
    \frac {\partial} {\partial t} \left[\begin{array}{c} \rho \\ \rho u \\ \rho v \\ \rho w \\ e \end{array}\right]
  + \frac {\partial} {\partial x} \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ \rho u w \\ (e+p) u\end{array}\right]
  + \frac {\partial} {\partial y} \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ \rho v w \\ (e+p) v \end{array}\right]
  + \frac {\partial} {\partial z} \left[\begin{array}{c} \rho w \\ \rho u w \\ \rho v w \\ \rho w^2 + p \\ (e+p) w \end{array}\right]
  = \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ \tau_{zx} \\ u \tau_{xx} + v \tau_{yx} + w \tau_{zx} - q_x \end{array}\right]
  + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{zy} \\ u \tau_{xy} + v \tau_{yy} + w \tau_{zy} - q_y \end{array}\right]
  + \frac {\partial} {\partial z} \left[\begin{array}{c} 0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\ u \tau_{xz} + v \tau_{yz} + w \tau_{zz} - q_z \end{array}\right]
  + \left[\begin{array}{c} 0 \\ -\rho {\bf g}\cdot{\bf \hat{i}} \\ -\rho {\bf g}\cdot{\bf \hat{j}}  \\ -\rho {\bf g}\cdot{\bf \hat{k}} \\ -\rho u {\bf g}\cdot{\bf \hat{i}} - \rho v {\bf g}\cdot{\bf \hat{j}} - \rho w {\bf g}\cdot{\bf \hat{k}} \end{array}\right]
  \f}
  where \f${\bf g}\f$ is the gravitational force vector per unit mass, \f${\bf \hat{i}},{\bf \hat{j}}, {\bf \hat{k}}\f$ are the unit vectors along the x, y, and z, the viscous terms are given by
  \f{align}{
    \tau_{ij} &= \frac{\mu}{Re_\infty} \left[ \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) - \frac{2}{3}\frac{\partial u_k}{\partial x_k} \delta_{ij} \right], \\
    q_i &= - \frac{\mu}{\left(\gamma-1\right)Re_\infty Pr} \frac{\partial T}{\partial x_i}
  \f}
  with \f$\mu\f$ being the viscosity coefficient (computed using Sutherland's law), and the equation of state is
  \f{equation}{
    e = \frac {p} {\gamma-1} + \frac{1}{2} \rho \left(u^2 + v^2 + w^2\right)
  \f}
  References for the governing equations (as well as non-dimensional form):-
  + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
    Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
    number is based on speed of sound, instead of the freestream velocity).

  Reference for the well-balanced treatment of gravitational source term:
  + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source 
    Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889, 
    7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
    http://dx.doi.org/10.2514/6.2015-2889
  + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm 
    for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.

  Reference for the partitioning of the flux into its stiff (acoustic) and non-stiff (convective)
  components:
  + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows 
    with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
    38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.

*/

#include <basic.h>

/*! 3D Navier-Stokes equations */
#define _NAVIER_STOKES_3D_  "navierstokes3d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 3
/*! Number of vector components at each grid point */
#define _MODEL_NVARS_ 5

/* choices for upwinding schemes */
/*! Roe's upwinding scheme */
#define _ROE_   "roe"
/*! Characteristic-based Roe-fixed scheme */
#define _RF_    "rf-char"
/*! Characteristic-based local Lax-Friedrich scheme */
#define _LLF_   "llf-char"
/*! Rusanov's upwinding scheme */
#define _RUSANOV_   "rusanov"

/* grid directions */
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0
/*! dimension corresponding to the \a y spatial dimension */
#define _YDIR_ 1
/*! dimension corresponding to the \a z spatial dimension */
#define _ZDIR_ 2

/*! \def _NavierStokes3DGetFlowVar_
 Get the flow variables from the conserved solution vector.
 \f{equation}{
   {\bf u} = \left[\begin{array}{c} \rho \\ \rho u \\ \rho v \\ \rho w \\ e \end{array}\right]
 \f}
*/
#define _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p) \
  { \
    double gamma   = p->gamma, vsq; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    vz  = u[3] / rho; \
    e   = u[4]; \
    vsq  = vx*vx + vy*vy + vz*vz; \
    P   = (e - 0.5*rho*vsq) * (gamma-1.0); \
  }

/*! \def _NavierStokes3DSetFlux_
  Compute the flux vector, given the flow variables
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ \rho u w \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ \rho v w \\ (e+p)v \end{array}\right], \\
    dir = z, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho w \\ \rho u w \\ \rho v w \\ \rho w^2 + p \\ (e+p)w \end{array}\right]
  \f}
*/
#define _NavierStokes3DSetFlux_(f,rho,vx,vy,vz,e,P,dir) \
  { \
    if (dir == _XDIR_) { \
      f[0] = rho * vx; \
      f[1] = rho * vx * vx + P; \
      f[2] = rho * vx * vy; \
      f[3] = rho * vx * vz; \
      f[4] = (e + P) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = rho * vy; \
      f[1] = rho * vy * vx; \
      f[2] = rho * vy * vy + P; \
      f[3] = rho * vy * vz; \
      f[4] = (e + P) * vy; \
    } else if (dir == _ZDIR_) { \
      f[0] = rho * vz; \
      f[1] = rho * vz * vx; \
      f[2] = rho * vz * vy; \
      f[3] = rho * vz * vz + P; \
      f[4] = (e + P) * vz; \
    } \
  }

/*! \def _NavierStokes3DSetStiffFlux_
  Compute the stiff flux vector (comprising the acoustic modes only), given the flow variables
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{1}{\gamma}\rho u \\ \frac{1}{\gamma}\rho u^2 + p \\ \frac{1}{\gamma}\rho u v \\ \frac{1}{\gamma}\rho u w \\ (e+p)u - \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2+w^2\right)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{1}{\gamma}\rho v \\ \frac{1}{\gamma}\rho u v \\ \frac{1}{\gamma}\rho v^2 + p \\ \frac{1}{\gamma}\rho v w \\ (e+p)v - \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2+w^2\right)v \end{array}\right], \\
    dir = z, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{1}{\gamma}\rho w \\ \frac{1}{\gamma}\rho u w \\ \frac{1}{\gamma}\rho v w \\ \frac{1}{\gamma}\rho w^2 + p \\ (e+p)w - \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2+w^2\right)w \end{array}\right], \\
  \f}
  Reference:
  + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
    with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
    38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
#define _NavierStokes3DSetStiffFlux_(f,rho,vx,vy,vz,e,P,dir,gamma) \
  { \
    double gamma_inv = 1.0/gamma; \
    if (dir == _XDIR_) { \
      f[0] = gamma_inv * rho * vx; \
      f[1] = gamma_inv * rho * vx * vx + P; \
      f[2] = gamma_inv * rho * vx * vy; \
      f[3] = gamma_inv * rho * vx * vz; \
      f[4] = (e + P) * vx - 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy+vz*vz) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = gamma_inv * rho * vy; \
      f[1] = gamma_inv * rho * vy * vx; \
      f[2] = gamma_inv * rho * vy * vy + P; \
      f[3] = gamma_inv * rho * vy * vz; \
      f[4] = (e + P) * vy - 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy+vz*vz) * vy; \
    } else if (dir == _ZDIR_) { \
      f[0] = gamma_inv * rho * vz; \
      f[1] = gamma_inv * rho * vz * vx; \
      f[2] = gamma_inv * rho * vz * vy; \
      f[3] = gamma_inv * rho * vz * vz + P; \
      f[4] = (e + P) * vz - 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy+vz*vz) * vz; \
    } \
  }

/*! \def _NavierStokes3DSetNonStiffFlux_
  Compute the non-stiff flux vector (comprising the entropy modes only), given the flow variables
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{\gamma-1}{\gamma}\rho u \\ \frac{\gamma-1}{\gamma}\rho u^2 \\ \frac{\gamma-11}{\gamma}\rho u v \\ \frac{\gamma-1}{\gamma}\rho u w \\ \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2+w^2\right)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{\gamma-1}{\gamma}\rho v \\ \frac{\gamma-1}{\gamma}\rho u v \\ \frac{\gamma-1}{\gamma}\rho v^2 \\ \frac{\gamma-1}{\gamma}\rho v w  \\ \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2+w^2\right)v \end{array}\right], \\
    dir = z, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{\gamma-1}{\gamma}\rho w \\ \frac{\gamma-1}{\gamma}\rho u w \\ \frac{\gamma-1}{\gamma}\rho v w \\ \frac{\gamma-1}{\gamma}\rho w^2  \\ \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2+w^2\right)w \end{array}\right], \\
  \f}
  Reference:
  + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
    with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
    38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
#define _NavierStokes3DSetNonStiffFlux_(f,rho,vx,vy,vz,e,P,dir,gamma) \
  { \
    double gamma_inv = 1.0/gamma; \
    if (dir == _XDIR_) { \
      f[0] = (gamma-1.0) * gamma_inv * rho * vx; \
      f[1] = (gamma-1.0) * gamma_inv * rho * vx * vx; \
      f[2] = (gamma-1.0) * gamma_inv * rho * vx * vy; \
      f[3] = (gamma-1.0) * gamma_inv * rho * vx * vz; \
      f[4] = 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy+vz*vz) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = (gamma-1.0) * gamma_inv * rho * vy; \
      f[1] = (gamma-1.0) * gamma_inv * rho * vy * vx; \
      f[2] = (gamma-1.0) * gamma_inv * rho * vy * vy; \
      f[3] = (gamma-1.0) * gamma_inv * rho * vy * vz; \
      f[4] = 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy+vz*vz) * vy; \
    } else if (dir == _ZDIR_) { \
      f[0] = (gamma-1.0) * gamma_inv * rho * vz; \
      f[1] = (gamma-1.0) * gamma_inv * rho * vz * vx; \
      f[2] = (gamma-1.0) * gamma_inv * rho * vz * vy; \
      f[3] = (gamma-1.0) * gamma_inv * rho * vz * vz; \
      f[4] = 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy+vz*vz) * vz; \
    } \
  }

/*! \def _NavierStokes3DRoeAverage_
  Compute the Roe-average of two solutions.
*/
#define _NavierStokes3DRoeAverage_(uavg,uL,uR,p) \
  { \
    double  rho ,vx, vy, vz, e ,P ,H; \
    double  rhoL,vxL,vyL,vzL,eL,PL,HL,cLsq; \
    double  rhoR,vxR,vyR,vzR,eR,PR,HR,cRsq; \
    double  gamma = p->gamma; \
    _NavierStokes3DGetFlowVar_(uL,rhoL,vxL,vyL,vzL,eL,PL,p); \
    cLsq = gamma * PL/rhoL; \
    HL = 0.5*(vxL*vxL+vyL*vyL+vzL*vzL) + cLsq / (gamma-1.0); \
    _NavierStokes3DGetFlowVar_(uR,rhoR,vxR,vyR,vzR,eR,PR,p); \
    cRsq = gamma * PR/rhoR; \
    HR = 0.5*(vxR*vxR+vyR*vyR+vzR*vzR) + cRsq / (gamma-1.0); \
    double tL = sqrt(rhoL); \
    double tR = sqrt(rhoR); \
    rho = tL * tR; \
    vx  = (tL*vxL + tR*vxR) / (tL + tR); \
    vy  = (tL*vyL + tR*vyR) / (tL + tR); \
    vz  = (tL*vzL + tR*vzR) / (tL + tR); \
    H   = (tL*HL  + tR*HR ) / (tL + tR); \
    P = (H - 0.5* (vx*vx+vy*vy+vz*vz)) * (rho*(gamma-1.0))/gamma; \
    e   = P/(gamma-1.0) + 0.5*rho*(vx*vx+vy*vy+vz*vz); \
    uavg[0] = rho; \
    uavg[1] = rho*vx; \
    uavg[2] = rho*vy; \
    uavg[3] = rho*vz; \
    uavg[4] = e; \
  }

/*! \def _NavierStokes3DEigenvalues_
  Compute the eigenvalues, given a solution vector in terms of the conserved variables. The eigenvalues are returned
  as a matrix D whose diagonal values are the eigenvalues. Admittedly, this is inefficient. The matrix D is stored in
  a row-major format.

  The eigenvalues are \f$ {\bf v}\cdot\hat{\bf n}, {\bf v}\cdot\hat{\bf n}, {\bf v}\cdot\hat{\bf n}, {\bf v}\cdot\hat{\bf n} \pm c\f$
  where \f${\bf v} = u\hat{\bf i} + v\hat{\bf j} + w\hat{\bf k}\f$ is the velocity vector, \f$\hat{\bf n} = \hat{\bf i}, \hat{\bf j}, \hat{\bf k}\f$
  is the unit vector along the spatial dimension, and \f$c\f$ is the speed of sound.
  \b Note that the order of the eigenvalues (and therefore, the order of the eigenvectors in #_NavierStokes3DLeftEigenvectors_ and 
  #_NavierStokes3DRightEigenvectors_ has been chosen to avoid left eigen-matrices with zero on the diagonals.
*/
#define _NavierStokes3DEigenvalues_(u,D,p,dir) \
  { \
    double          gamma   = p->gamma; \
    double          rho,vx,vy,vz,e,P,c; \
    _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p); \
    c    = sqrt(gamma*P/rho); \
    if (dir == _XDIR_) { \
      D[0*_MODEL_NVARS_+0] = vx;     D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0;    D[0*_MODEL_NVARS_+4] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vx-c; D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0;    D[1*_MODEL_NVARS_+4] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vx;     D[2*_MODEL_NVARS_+3] = 0;    D[2*_MODEL_NVARS_+4] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vx;   D[3*_MODEL_NVARS_+4] = 0; \
      D[4*_MODEL_NVARS_+0] = 0;      D[4*_MODEL_NVARS_+1] = 0;    D[4*_MODEL_NVARS_+2] = 0;      D[4*_MODEL_NVARS_+3] = 0;    D[4*_MODEL_NVARS_+4] = vx+c;\
    } else if (dir == _YDIR_) { \
      D[0*_MODEL_NVARS_+0] = vy;     D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0;    D[0*_MODEL_NVARS_+4] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vy;   D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0;    D[1*_MODEL_NVARS_+4] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vy-c;   D[2*_MODEL_NVARS_+3] = 0;    D[2*_MODEL_NVARS_+4] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vy;   D[3*_MODEL_NVARS_+4] = 0; \
      D[4*_MODEL_NVARS_+0] = 0;      D[4*_MODEL_NVARS_+1] = 0;    D[4*_MODEL_NVARS_+2] = 0;      D[4*_MODEL_NVARS_+3] = 0;    D[4*_MODEL_NVARS_+4] = vy+c;\
    } else if (dir == _ZDIR_) { \
      D[0*_MODEL_NVARS_+0] = vz;     D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0;    D[0*_MODEL_NVARS_+4] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vz;   D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0;    D[1*_MODEL_NVARS_+4] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vz;     D[2*_MODEL_NVARS_+3] = 0;    D[2*_MODEL_NVARS_+4] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vz-c; D[3*_MODEL_NVARS_+4] = 0; \
      D[4*_MODEL_NVARS_+0] = 0;      D[4*_MODEL_NVARS_+1] = 0;    D[4*_MODEL_NVARS_+2] = 0;      D[4*_MODEL_NVARS_+3] = 0;    D[4*_MODEL_NVARS_+4] = vz+c;\
    } \
  }

/*! \def _NavierStokes3DLeftEigenvectors_
  Compute the left eigenvectors, given a solution vector in terms of the conserved variables. The eigenvectors are
  returned as a matrix L whose rows correspond to each eigenvector. The matrix L is stored in the row-major format.
  \n\n
  Reference:
  + Rohde, "Eigenvalues and eigenvectors of the Euler equations in general geometries", AIAA Paper 2001-2609,
    http://dx.doi.org/10.2514/6.2001-2609
*/
#define _NavierStokes3DLeftEigenvectors_(u,L,p,dir) \
  { \
    double  ga = p->gamma, ga_minus_one=ga-1.0; \
    double  rho,vx,vy,vz,e,P,a,ek; \
    _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p); \
  	ek = 0.5 * (vx*vx + vy*vy + vz*vz); \
	  a = sqrt(ga * P / rho); \
    if (dir == _XDIR_) { \
    	L[1*_MODEL_NVARS_+0] = (ga_minus_one*ek + a*vx) / (2*a*a); \
  		L[1*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx-a) / (2*a*a); \
  		L[1*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy) / (2*a*a); \
  		L[1*_MODEL_NVARS_+3] = ((-ga_minus_one)*vz) / (2*a*a); \
  		L[1*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[0*_MODEL_NVARS_+0] = (a*a - ga_minus_one*ek) / (a*a); \
  		L[0*_MODEL_NVARS_+1] = (ga_minus_one*vx) / (a*a); \
  		L[0*_MODEL_NVARS_+2] = (ga_minus_one*vy) / (a*a); \
  		L[0*_MODEL_NVARS_+3] = (ga_minus_one*vz) / (a*a); \
  		L[0*_MODEL_NVARS_+4] = (-ga_minus_one) / (a*a); \
  		L[4*_MODEL_NVARS_+0] = (ga_minus_one*ek - a*vx) / (2*a*a); \
  		L[4*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx+a) / (2*a*a); \
  		L[4*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy) / (2*a*a); \
  		L[4*_MODEL_NVARS_+3] = ((-ga_minus_one)*vz) / (2*a*a); \
  		L[4*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[2*_MODEL_NVARS_+0] = vy; \
  		L[2*_MODEL_NVARS_+1] = 0.0; \
  		L[2*_MODEL_NVARS_+2] = -1.0; \
  		L[2*_MODEL_NVARS_+3] = 0.0; \
  		L[2*_MODEL_NVARS_+4] = 0.0; \
  		L[3*_MODEL_NVARS_+0] = -vz; \
  		L[3*_MODEL_NVARS_+1] = 0.0; \
  		L[3*_MODEL_NVARS_+2] = 0.0; \
  		L[3*_MODEL_NVARS_+3] = 1.0; \
  		L[3*_MODEL_NVARS_+4] = 0.0; \
    } else if (dir == _YDIR_) {  \
  		L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*vy) / (2*a*a); \
  		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx) / (2*a*a); \
  		L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy-a) / (2*a*a); \
  		L[2*_MODEL_NVARS_+3] = ((1.0-ga)*vz) / (2*a*a); \
  		L[2*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[0*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a); \
  		L[0*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a); \
  		L[0*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a); \
  		L[0*_MODEL_NVARS_+3] = ga_minus_one*vz / (a*a); \
  		L[0*_MODEL_NVARS_+4] = (1.0 - ga) / (a*a); \
  		L[4*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*vy) / (2*a*a); \
  		L[4*_MODEL_NVARS_+1] = ((1.0-ga)*vx) / (2*a*a); \
  		L[4*_MODEL_NVARS_+2] = ((1.0-ga)*vy+a) / (2*a*a); \
  		L[4*_MODEL_NVARS_+3] = ((1.0-ga)*vz) / (2*a*a); \
  		L[4*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[1*_MODEL_NVARS_+0] = -vx; \
  		L[1*_MODEL_NVARS_+1] = 1.0; \
  		L[1*_MODEL_NVARS_+2] = 0.0; \
  		L[1*_MODEL_NVARS_+3] = 0.0; \
  		L[1*_MODEL_NVARS_+4] = 0; \
  		L[3*_MODEL_NVARS_+0] = vz; \
  		L[3*_MODEL_NVARS_+1] = 0.0; \
  		L[3*_MODEL_NVARS_+2] = 0.0; \
  		L[3*_MODEL_NVARS_+3] = -1.0; \
  		L[3*_MODEL_NVARS_+4] = 0; \
    } else if (dir == _ZDIR_) {  \
      L[3*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*vz) / (2*a*a); \
  		L[3*_MODEL_NVARS_+1] = ((1.0-ga)*vx) / (2*a*a); \
  		L[3*_MODEL_NVARS_+2] = ((1.0-ga)*vy) / (2*a*a); \
  		L[3*_MODEL_NVARS_+3] = ((1.0-ga)*vz-a) / (2*a*a); \
  		L[3*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[0*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a); \
  		L[0*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a); \
  		L[0*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a); \
  		L[0*_MODEL_NVARS_+3] = ga_minus_one*vz / (a*a); \
  		L[0*_MODEL_NVARS_+4] = (1.0-ga) / (a*a); \
  		L[4*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*vz) / (2*a*a); \
  		L[4*_MODEL_NVARS_+1] = ((1.0-ga)*vx) / (2*a*a); \
  		L[4*_MODEL_NVARS_+2] = ((1.0-ga)*vy) / (2*a*a); \
  		L[4*_MODEL_NVARS_+3] = ((1.0-ga)*vz+a) / (2*a*a); \
  		L[4*_MODEL_NVARS_+4] = ga_minus_one / (2*a*a); \
  		L[1*_MODEL_NVARS_+0] = vx; \
  		L[1*_MODEL_NVARS_+1] = -1.0; \
  		L[1*_MODEL_NVARS_+2] = 0.0; \
  		L[1*_MODEL_NVARS_+3] = 0.0; \
  		L[1*_MODEL_NVARS_+4] = 0; \
  		L[2*_MODEL_NVARS_+0] = -vy; \
  		L[2*_MODEL_NVARS_+1] = 0.0; \
  		L[2*_MODEL_NVARS_+2] = 1.0; \
  		L[2*_MODEL_NVARS_+3] = 0.0; \
  		L[2*_MODEL_NVARS_+4] = 0; \
    } \
  }

/*! \def _NavierStokes3DRightEigenvectors_
  Compute the right eigenvectors, given a solution vector in terms of the conserved variables. The eigenvectors are
  returned as a matrix R whose columns correspond to each eigenvector. The matrix R is stored in the row-major format.
  \n\n
  Reference:
  + Rohde, "Eigenvalues and eigenvectors of the Euler equations in general geometries", AIAA Paper 2001-2609,
    http://dx.doi.org/10.2514/6.2001-2609
*/
#define _NavierStokes3DRightEigenvectors_(u,R,p,dir) \
  { \
    double  ga   = p->gamma, ga_minus_one = ga-1.0; \
    double  rho,vx,vy,vz,e,P,ek,a,h0; \
    _NavierStokes3DGetFlowVar_(u,rho,vx,vy,vz,e,P,p); \
  	ek   = 0.5 * (vx*vx + vy*vy + vz*vz); \
	  a    = sqrt(ga * P / rho); \
    h0   = a*a / ga_minus_one + ek; \
	  if (dir == _XDIR_) { \
	  	R[0*_MODEL_NVARS_+1] = 1.0; \
		  R[1*_MODEL_NVARS_+1] = vx-a; \
  		R[2*_MODEL_NVARS_+1] = vy; \
	  	R[3*_MODEL_NVARS_+1] = vz; \
		  R[4*_MODEL_NVARS_+1] = h0 - a*vx; \
  		R[0*_MODEL_NVARS_+0] = 1.0; \
	  	R[1*_MODEL_NVARS_+0] = vx; \
		  R[2*_MODEL_NVARS_+0] = vy; \
  		R[3*_MODEL_NVARS_+0] = vz; \
	  	R[4*_MODEL_NVARS_+0] = ek; \
		  R[0*_MODEL_NVARS_+4] = 1.0; \
  		R[1*_MODEL_NVARS_+4] = vx+a; \
	  	R[2*_MODEL_NVARS_+4] = vy; \
		  R[3*_MODEL_NVARS_+4] = vz; \
  		R[4*_MODEL_NVARS_+4] = h0 + a*vx; \
	  	R[0*_MODEL_NVARS_+2] = 0.0; \
		  R[1*_MODEL_NVARS_+2] = 0.0; \
  		R[2*_MODEL_NVARS_+2] = -1.0; \
	  	R[3*_MODEL_NVARS_+2] = 0.0; \
		  R[4*_MODEL_NVARS_+2] = -vy; \
  		R[0*_MODEL_NVARS_+3] = 0.0; \
	  	R[1*_MODEL_NVARS_+3] = 0.0; \
		  R[2*_MODEL_NVARS_+3] = 0.0; \
  		R[3*_MODEL_NVARS_+3] = 1.0; \
	  	R[4*_MODEL_NVARS_+3] = vz; \
	  } else if (dir == _YDIR_) { \
	  	R[0*_MODEL_NVARS_+2] = 1.0; \
  		R[1*_MODEL_NVARS_+2] = vx; \
	  	R[2*_MODEL_NVARS_+2] = vy-a; \
		  R[3*_MODEL_NVARS_+2] = vz; \
  		R[4*_MODEL_NVARS_+2] = h0 - a*vy; \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
  		R[1*_MODEL_NVARS_+0] = vx; \
	  	R[2*_MODEL_NVARS_+0] = vy; \
		  R[3*_MODEL_NVARS_+0] = vz; \
  		R[4*_MODEL_NVARS_+0] = ek; \
	  	R[0*_MODEL_NVARS_+4] = 1.0; \
		  R[1*_MODEL_NVARS_+4] = vx; \
  		R[2*_MODEL_NVARS_+4] = vy+a; \
	  	R[3*_MODEL_NVARS_+4] = vz; \
		  R[4*_MODEL_NVARS_+4] = h0 + a*vy; \
  		R[0*_MODEL_NVARS_+1] = 0.0; \
	  	R[1*_MODEL_NVARS_+1] = 1.0; \
		  R[2*_MODEL_NVARS_+1] = 0.0; \
  		R[3*_MODEL_NVARS_+1] = 0.0; \
	  	R[4*_MODEL_NVARS_+1] = vx; \
  		R[0*_MODEL_NVARS_+3] = 0.0; \
	  	R[1*_MODEL_NVARS_+3] = 0.0; \
		  R[2*_MODEL_NVARS_+3] = 0.0; \
  		R[3*_MODEL_NVARS_+3] = -1.0; \
	  	R[4*_MODEL_NVARS_+3] = -vz; \
    } else if (dir == _ZDIR_) {  \
	  	R[0*_MODEL_NVARS_+3] = 1.0; \
	  	R[1*_MODEL_NVARS_+3] = vx; \
	  	R[2*_MODEL_NVARS_+3] = vy; \
	  	R[3*_MODEL_NVARS_+3] = vz-a; \
	  	R[4*_MODEL_NVARS_+3] = h0-a*vz; \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
	  	R[1*_MODEL_NVARS_+0] = vx; \
	  	R[2*_MODEL_NVARS_+0] = vy; \
	  	R[3*_MODEL_NVARS_+0] = vz; \
	  	R[4*_MODEL_NVARS_+0] = ek; \
	  	R[0*_MODEL_NVARS_+4] = 1.0; \
	  	R[1*_MODEL_NVARS_+4] = vx; \
	  	R[2*_MODEL_NVARS_+4] = vy; \
	  	R[3*_MODEL_NVARS_+4] = vz+a; \
	  	R[4*_MODEL_NVARS_+4] = h0+a*vz; \
	  	R[0*_MODEL_NVARS_+1] = 0.0; \
	  	R[1*_MODEL_NVARS_+1] = -1.0; \
	  	R[2*_MODEL_NVARS_+1] = 0.0; \
	  	R[3*_MODEL_NVARS_+1] = 0.0; \
	  	R[4*_MODEL_NVARS_+1] = -vx; \
	  	R[0*_MODEL_NVARS_+2] = 0.0; \
	  	R[1*_MODEL_NVARS_+2] = 0.0; \
	  	R[2*_MODEL_NVARS_+2] = 1.0; \
	  	R[3*_MODEL_NVARS_+2] = 0.0; \
	  	R[4*_MODEL_NVARS_+2] = vy; \
    } \
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
  double  gamma;  /*!< Ratio of heat capacities */
  double  Re;     /*!< Reynolds number */
  double  Pr;     /*!< Prandtl  number */
  double  Minf;   /*!< Freestream Mach number */
  double  C1,     /*!< Sutherlands law constant */
          C2;     /*!< Sutherlands law constant */
  double  grav_x,                         /*!< acceleration due to gravity in x */
          grav_y,                         /*!< acceleration due to gravity in y */
          grav_z;                         /*!< acceleration due to gravity in z */
  double  rho0;                           /*!< reference density  at zero altitude for flows with gravity */
  double  p0;                             /*!< reference pressure at zero altitude for flows with gravity */
  double  R;                              /*!< universal Gas constant */
  char    upw_choice[_MAX_STRING_SIZE_]; /*!< choice of upwinding */

  double  *grav_field_f, /*!< density variation function (\f$\varrho\f$) for hydrostatic equilibrium for flows with gravity */
          *grav_field_g; /*!< pressure variation function (\f$\varrho\f$) for hydrostatic equilibrium for flows with gravity */

  double *fast_jac, /*!< "Fast" Jacobian of the flux function (comprising the acoustic modes) */
         *solution; /*!< array to store the solution at the beginning of each time step */

  /* choice of hydrostatic balance                              */
  /* 1 -> isothermal                                            */
  /* 2 -> constant potential temperature                        */ 
  /* 3 -> stratified atmosphere with a Brunt-Vaisala frequency  */
  int HB /*!< Choice of hydrostatic balance for flows with gravity (1 - isothermal equilibrium, 
                                                                    2 - constant potential temperature
                                                                    3 - stratified atmosphere with a Brunt-Vaisala frequency) */;
  double N_bv; /*!< the Brunt-Vaisala frequency for #NavierStokes3D::HB = 3 */

  char ib_write_surface_data[_MAX_STRING_SIZE_]; /*!< Flag to indicate whether to analyze and write surface data for
                                                      immersed body, if present. Applicable only if #HyPar::flag_ib is 1 */

} NavierStokes3D;

int    NavierStokes3DInitialize (void*,void*);
int    NavierStokes3DCleanup    (void*);

