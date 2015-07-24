/*! @file navierstokes2d.h
    @brief 2D Navier Stokes equations (compressible flows)
    @author Debojyoti Ghosh

  2D Navier-Stokes equations for viscous and inviscid compressible flows (with gravitational terms)\n

  \f{equation}{
    \frac {\partial} {\partial t} \left[\begin{array}{c} \rho \\ \rho u \\ \rho v \\ e \end{array}\right]
  + \frac {\partial} {\partial x} \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ (e+p) u\end{array}\right]
  + \frac {\partial} {\partial y} \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ (e+p) v \end{array}\right]
  = \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ u \tau_{xx} + v \tau_{yx} - q_x \end{array}\right]
  + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ u \tau_{xy} + v \tau_{yy} - q_y \end{array}\right]
  + \left[\begin{array}{c} 0 \\ -\rho {\bf g}\cdot{\bf \hat{i}} \\ -\rho {\bf g}\cdot{\bf \hat{j}}  \\ -\rho u {\bf g}\cdot{\bf \hat{i}} - \rho v {\bf g}\cdot{\bf \hat{j}} \end{array}\right]
  \f}
  where \f${\bf g}\f$ is the gravitational force vector per unit mass, \f${\bf \hat{i}},{\bf \hat{j}}\f$ are the unit vectors along the x and y, the viscous terms are given by
  \f{align}{
    \tau_{ij} &= \frac{\mu}{Re_\infty} \left[ \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) - \frac{2}{3}\frac{\partial u_k}{\partial x_k} \delta_{ij} \right], \\
    q_i &= - \frac{\mu}{\left(\gamma-1\right)Re_\infty Pr} \frac{\partial T}{\partial x_i}
  \f}
  with \f$\mu\f$ being the viscosity coefficient (computed using Sutherland's law), and the equation of state is
  \f{equation}{
    e = \frac {p} {\gamma-1} + \frac{1}{2} \rho \left(u^2 + v^2\right)
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
    for Atmospheric Flows, Submitted
*/
#include <basic.h>

/*! 2D Navier Stokes equations */
#define _NAVIER_STOKES_2D_  "navierstokes2d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 2
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 4

/* choices for upwinding schemes */
/*! Roe's upwinding scheme */
#define _ROE_       "roe"
/*! Characteristic-based Roe-fixed scheme */
#define _RF_        "rf-char"
/*! Characteristic-based local Lax-Friedrich scheme */
#define _LLF_       "llf-char"
/*! Steger-Warming flux splitting scheme */
#define _SWFS_      "steger-warming"
/*! Rusanov's upwinding scheme */
#define _RUSANOV_   "rusanov"

/* directions */
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0
/*! dimension corresponding to the \a y spatial dimension */
#define _YDIR_ 1

/*! \def _NavierStokes2DGetFlowVar_
 Get the flow variables from the conserved solution vector.
 \f{equation}{
   {\bf u} = \left[\begin{array}{c} \rho \\ \rho u \\ \rho v \\ e \end{array}\right]
 \f}
*/
#define _NavierStokes2DGetFlowVar_(u,rho,vx,vy,e,P,p) \
  { \
    double  gamma = p->gamma, vsq; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (gamma-1.0); \
  }

/*! \def _NavierStokes2DSetFlux_
  Compute the flux vector, given the flow variables
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ (e+p)v \end{array}\right]
  \f}
*/
#define _NavierStokes2DSetFlux_(f,rho,vx,vy,e,P,dir) \
  { \
    if (dir == _XDIR_) { \
      f[0] = rho * vx; \
      f[1] = rho * vx * vx + P; \
      f[2] = rho * vx * vy; \
      f[3] = (e + P) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = rho * vy; \
      f[1] = rho * vy * vx; \
      f[2] = rho * vy * vy + P; \
      f[3] = (e + P) * vy; \
    } \
  }

/*! \def _NavierStokes2DSetStiffFlux_
  Compute the stiff flux vector (comprising the acoustic modes only), given the flow variables
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{1}{\gamma}\rho u \\ \frac{1}{\gamma}\rho u^2 + p \\ \frac{1}{\gamma}\rho u v \\ (e+p)u - \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2\right)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \frac{1}{\gamma}\rho v \\ \frac{1}{\gamma}\rho u v \\ \frac{1}{\gamma}\rho v^2 + p \\ (e+p)v - \frac{1}{2} \frac{\gamma-1}{\gamma}\rho\left(u^2+v^2\right)v \end{array}\right]
  \f}
  Reference:
  + Not yet published!
*/
#define _NavierStokes2DSetStiffFlux_(f,rho,vx,vy,e,P,dir,gamma) \
  { \
    double gamma_inv = 1.0/gamma; \
    if (dir == _XDIR_) { \
      f[0] = gamma_inv * rho * vx; \
      f[1] = gamma_inv * rho * vx * vx + P; \
      f[2] = gamma_inv * rho * vx * vy; \
      f[3] = (e + P) * vx - 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy) * vx; \
    } else if (dir == _YDIR_) { \
      f[0] = gamma_inv * rho * vy; \
      f[1] = gamma_inv * rho * vy * vx; \
      f[2] = gamma_inv * rho * vy * vy + P; \
      f[3] = (e + P) * vy - 0.5 * gamma_inv * (gamma-1.0) * rho * (vx*vx+vy*vy) * vy; \
    } \
  }

/*! \def _NavierStokes2DRoeAverage_
  Compute the Roe-average of two solutions.
*/
#define _NavierStokes2DRoeAverage_(uavg,uL,uR,p) \
  { \
    double  rho ,vx, vy, e ,P ,H ,csq, vsq; \
    double  rhoL,vxL,vyL,eL,PL,HL,cLsq,vsqL; \
    double  rhoR,vxR,vyR,eR,PR,HR,cRsq,vsqR; \
    double  gamma = p->gamma; \
    rhoL = uL[0]; \
    vxL  = uL[1] / rhoL; \
    vyL  = uL[2] / rhoL; \
    eL   = uL[3]; \
    vsqL = (vxL*vxL) + (vyL*vyL); \
    PL   = (eL - 0.5*rhoL*vsqL) * (gamma-1.0); \
    cLsq = gamma * PL/rhoL; \
    HL = 0.5*(vxL*vxL+vyL*vyL) + cLsq / (gamma-1.0); \
    rhoR = uR[0]; \
    vxR  = uR[1] / rhoR; \
    vyR  = uR[2] / rhoR; \
    eR   = uR[3]; \
    vsqR = (vxR*vxR) + (vyR*vyR); \
    PR   = (eR - 0.5*rhoR*vsqR) * (gamma-1.0); \
    cRsq = gamma * PR/rhoR; \
    HR = 0.5*(vxR*vxR+vyR*vyR) + cRsq / (gamma-1.0); \
    double tL = sqrt(rhoL); \
    double tR = sqrt(rhoR); \
    rho = tL * tR; \
    vx  = (tL*vxL + tR*vxR) / (tL + tR); \
    vy  = (tL*vyL + tR*vyR) / (tL + tR); \
    H   = (tL*HL + tR*HR) / (tL + tR); \
    vsq = vx*vx + vy*vy; \
    csq = (gamma-1.0) * (H-0.5*vsq); \
    P   = csq * rho / gamma; \
    e   = P/(gamma-1.0) + 0.5*rho*vsq; \
    uavg[0] = rho; \
    uavg[1] = rho*vx; \
    uavg[2] = rho*vy; \
    uavg[3] = e; \
  }

/*! \def _NavierStokes2DEigenvalues_
  Compute the eigenvalues, given a solution vector in terms of the conserved variables. The eigenvalues are returned
  as a matrix D whose diagonal values are the eigenvalues. Admittedly, this is inefficient. The matrix D is stored in
  a row-major format.
*/
#define _NavierStokes2DEigenvalues_(u,D,p,dir) \
  { \
    double  gamma = p->gamma; \
    double  rho,vx,vy,e,P,c,vn,vsq; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (gamma-1.0); \
    c    = sqrt(gamma*P/rho); \
    if      (dir == _XDIR_) vn = vx; \
    else if (dir == _YDIR_) vn = vy; \
    else               vn = 0.0; \
    if (dir == _XDIR_) {\
      D[0*_MODEL_NVARS_+0] = vn-c;   D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vn+c; D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vn;     D[2*_MODEL_NVARS_+3] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vn; \
    } else if (dir == _YDIR_) { \
      D[0*_MODEL_NVARS_+0] = vn-c;   D[0*_MODEL_NVARS_+1] = 0;    D[0*_MODEL_NVARS_+2] = 0;      D[0*_MODEL_NVARS_+3] = 0; \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vn;   D[1*_MODEL_NVARS_+2] = 0;      D[1*_MODEL_NVARS_+3] = 0; \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;    D[2*_MODEL_NVARS_+2] = vn+c;   D[2*_MODEL_NVARS_+3] = 0; \
      D[3*_MODEL_NVARS_+0] = 0;      D[3*_MODEL_NVARS_+1] = 0;    D[3*_MODEL_NVARS_+2] = 0;      D[3*_MODEL_NVARS_+3] = vn; \
    }\
  }

/*! \def _NavierStokes2DLeftEigenvectors_
  Compute the left eigenvectors, given a solution vector in terms of the conserved variables. The eigenvectors are
  returned as a matrix L whose rows correspond to each eigenvector. The matrix L is stored in the row-major format.
  \n\n
  Reference:
  + Rohde, "Eigenvalues and eigenvectors of the Euler equations in general geometries", AIAA Paper 2001-2609,
    http://dx.doi.org/10.2514/6.2001-2609
*/
#define _NavierStokes2DLeftEigenvectors_(u,L,p,dir) \
  { \
    double  ga = param->gamma, ga_minus_one=ga-1.0; \
    double  rho,vx,vy,e,P,a,un,ek,vsq; \
    double  nx = 0,ny = 0; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (ga-1.0); \
  	ek = 0.5 * (vx*vx + vy*vy); \
	  a = sqrt(ga * P / rho); \
    if (dir == _XDIR_) { \
      un = vx; \
      nx = 1.0; \
  		L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek + a*un) / (2*a*a); \
	  	L[0*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx - a*nx) / (2*a*a); \
		  L[0*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy - a*ny) / (2*a*a); \
  		L[0*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
	  	L[3*_MODEL_NVARS_+0] = (a*a - ga_minus_one*ek) / (a*a); \
  		L[3*_MODEL_NVARS_+1] = (ga_minus_one*vx) / (a*a); \
	  	L[3*_MODEL_NVARS_+2] = (ga_minus_one*vy) / (a*a); \
		  L[3*_MODEL_NVARS_+3] = (-ga_minus_one) / (a*a); \
		  L[1*_MODEL_NVARS_+0] = (ga_minus_one*ek - a*un) / (2*a*a); \
  		L[1*_MODEL_NVARS_+1] = ((-ga_minus_one)*vx + a*nx) / (2*a*a); \
	  	L[1*_MODEL_NVARS_+2] = ((-ga_minus_one)*vy + a*ny) / (2*a*a); \
		  L[1*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
  		L[2*_MODEL_NVARS_+0] = (vy - un*ny) / nx; \
	  	L[2*_MODEL_NVARS_+1] = ny; \
		  L[2*_MODEL_NVARS_+2] = (ny*ny - 1.0) / nx; \
  		L[2*_MODEL_NVARS_+3] = 0.0; \
    } else if (dir == _YDIR_) {  \
      un = vy;  \
      ny = 1.0; \
	  	L[0*_MODEL_NVARS_+0] = (ga_minus_one*ek+a*un) / (2*a*a); \
		  L[0*_MODEL_NVARS_+1] = ((1.0-ga)*vx - a*nx) / (2*a*a); \
  		L[0*_MODEL_NVARS_+2] = ((1.0-ga)*vy - a*ny) / (2*a*a); \
	  	L[0*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
		  L[3*_MODEL_NVARS_+0] = (a*a-ga_minus_one*ek) / (a*a); \
  		L[3*_MODEL_NVARS_+1] = ga_minus_one*vx / (a*a); \
	  	L[3*_MODEL_NVARS_+2] = ga_minus_one*vy / (a*a); \
		  L[3*_MODEL_NVARS_+3] = (1.0 - ga) / (a*a); \
		  L[2*_MODEL_NVARS_+0] = (ga_minus_one*ek-a*un) / (2*a*a); \
  		L[2*_MODEL_NVARS_+1] = ((1.0-ga)*vx + a*nx) / (2*a*a); \
	  	L[2*_MODEL_NVARS_+2] = ((1.0-ga)*vy + a*ny) / (2*a*a); \
		  L[2*_MODEL_NVARS_+3] = ga_minus_one / (2*a*a); \
		  L[1*_MODEL_NVARS_+0] = (un*nx-vx) / ny; \
  		L[1*_MODEL_NVARS_+1] = (1.0 - nx*nx) / ny; \
	  	L[1*_MODEL_NVARS_+2] = - nx; \
		  L[1*_MODEL_NVARS_+3] = 0; \
    } \
  }

/*! \def _NavierStokes2DRightEigenvectors_
  Compute the right eigenvectors, given a solution vector in terms of the conserved variables. The eigenvectors are
  returned as a matrix R whose columns correspond to each eigenvector. The matrix R is stored in the row-major format.
  \n\n
  Reference:
  + Rohde, "Eigenvalues and eigenvectors of the Euler equations in general geometries", AIAA Paper 2001-2609,
    http://dx.doi.org/10.2514/6.2001-2609
*/
#define _NavierStokes2DRightEigenvectors_(u,R,p,dir) \
  { \
    double  ga   = param->gamma, ga_minus_one = ga-1.0; \
    double  rho,vx,vy,e,P,un,ek,a,h0,vsq; \
    double  nx = 0,ny = 0; \
    rho = u[0]; \
    vx  = u[1] / rho; \
    vy  = u[2] / rho; \
    e   = u[3]; \
    vsq  = (vx*vx) + (vy*vy); \
    P   = (e - 0.5*rho*vsq) * (ga-1.0); \
	  ek   = 0.5 * (vx*vx + vy*vy); \
  	a    = sqrt(ga * P / rho); \
    h0   = a*a / ga_minus_one + ek; \
	  if (dir == _XDIR_) { \
    	un = vx; \
      nx = 1.0; \
	  	R[0*_MODEL_NVARS_+0] = 1.0; \
	  	R[1*_MODEL_NVARS_+0] = vx - a*nx; \
	  	R[2*_MODEL_NVARS_+0] = vy - a*ny; \
	  	R[3*_MODEL_NVARS_+0] = h0 - a*un; \
  		R[0*_MODEL_NVARS_+3] = 1.0; \
  		R[1*_MODEL_NVARS_+3] = vx; \
  		R[2*_MODEL_NVARS_+3] = vy; \
  		R[3*_MODEL_NVARS_+3] = ek; \
  		R[0*_MODEL_NVARS_+1] = 1.0; \
  		R[1*_MODEL_NVARS_+1] = vx + a*nx; \
  		R[2*_MODEL_NVARS_+1] = vy + a*ny; \
  		R[3*_MODEL_NVARS_+1] = h0 + a*un; \
  		R[0*_MODEL_NVARS_+2] = 0.0; \
  		R[1*_MODEL_NVARS_+2] = ny; \
  		R[2*_MODEL_NVARS_+2] = -nx; \
  		R[3*_MODEL_NVARS_+2] = vx*ny - vy*nx; \
  	} else if (dir == _YDIR_) { \
      un = vy; \
      ny = 1.0; \
  		R[0*_MODEL_NVARS_+0] = 1.0; \
  		R[1*_MODEL_NVARS_+0] = vx - a*nx; \
  		R[2*_MODEL_NVARS_+0] = vy - a*ny; \
  		R[3*_MODEL_NVARS_+0] = h0 - a*un; \
  		R[0*_MODEL_NVARS_+3] = 1.0; \
  		R[1*_MODEL_NVARS_+3] = vx; \
  		R[2*_MODEL_NVARS_+3] = vy; \
  		R[3*_MODEL_NVARS_+3] = ek; \
  		R[0*_MODEL_NVARS_+2] = 1.0; \
  		R[1*_MODEL_NVARS_+2] = vx + a*nx; \
  		R[2*_MODEL_NVARS_+2] = vy + a*ny; \
  		R[3*_MODEL_NVARS_+2] = h0 + a*un; \
  		R[0*_MODEL_NVARS_+1] = 0; \
  		R[1*_MODEL_NVARS_+1] = ny; \
  		R[2*_MODEL_NVARS_+1] = -nx; \
  		R[3*_MODEL_NVARS_+1] = vx*ny-vy*nx; \
    } \
  }


/*! \def NavierStokes2D
    \brief Structure containing variables and parameters specific to the 2D Navier Stokes equations.
 *  This structure contains the physical parameters, variables, and function pointers specific to 
 *  the 2D Navier-Stokes equations.
*/
/*! \brief Structure containing variables and parameters specific to the 2D Navier Stokes equations.
 *  This structure contains the physical parameters, variables, and function pointers specific to 
 *  the 2D Navier-Stokes equations.
*/
typedef struct navierstokes2d_parameters {
  double  gamma;                          /*!< Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_];  /*!< choice of upwinding */
  double  grav_x,                         /*!< acceleration due to gravity in x */
          grav_y;                         /*!< acceleration due to gravity in y */
  double  rho0;                           /*!< reference density  at zero altitude for flows with gravity */
  double  p0;                             /*!< reference pressure at zero altitude for flows with gravity */
  double  Re;                             /*!< Reynolds number */
  double  Pr;                             /*!< Prandtl  number */
  double  Minf;                           /*!< Freestream Mach number */
  double  C1,                             /*!< Sutherlands law constants */
          C2;                             /*!< Sutherlands law constants */
  double  R;                              /*!< universal Gas constant */
  
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
  double N_bv; /*!< the Brunt-Vaisala frequency for #NavierStokes2D::HB = 3 */

} NavierStokes2D;

int    NavierStokes2DInitialize (void*,void*);
int    NavierStokes2DCleanup    (void*);

