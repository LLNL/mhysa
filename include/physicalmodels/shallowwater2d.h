/*! @file shallowwater2d.h
    @brief 2D Shallow Water Equations
    @author Debojyoti Ghosh

    2D Shallow Water Equations\n

    \f{equation}{
      \frac {\partial} {\partial t} \left[\begin{array}{c} h \\ hu \\ hv \end{array} \right]
    + \frac {\partial} {\partial x} \left[\begin{array}{c} hu \\ hu^2 + \frac{1}{2}gh^2 \\ huv  \end{array} \right] 
    + \frac {\partial} {\partial x} \left[\begin{array}{c} hv \\ huv \\ hv^2 + \frac{1}{2}gh^2  \end{array} \right] 
    = \left[\begin{array}{c} 0 \\ -ghb_x \\ -ghb_y \end{array}\right]
    + \left[\begin{array}{c} 0 \\ fv     \\ -fu    \end{array}\right]
    \f}
    where \f$h\f$ is the height, \f$\left(u,v\right)\f$ are the velocity components, \f$b(x,y)\f$ is the bottom 
    topography, and \f$g\f$ is the gravitational constant.\n\n

    The Coriolis parameter \f$f\f$ is defined as
    \f{equation}{
      f = \hat{f} + \beta\left( y - \frac{D}{2} \right)
    \f}
    based on:
    + Zhu, Et. al., "Variational Data Assimilation with a Variable Resolution
      Finite-Element Shallow-Water Equations Model", Monthly Weather Review,
      122, 1994, pp. 946--965
      http://dx.doi.org/10.1175/1520-0493(1994)122%3C0946:VDAWAV%3E2.0.CO;2\n
      Eqns. (2.1)-(2.3), (2.4)

    For the treatment of the topography gradient source term (well-balanced formulation), 
    refer to:\n
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the 
      exact conservation property for the shallow water equations", Journal of 
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/

#include <basic.h>
#include <math.h>
#include <matops.h>

/*! \def _SHALLOW_WATER_2D_
    2D Shallow Water equations
*/
#define _SHALLOW_WATER_2D_ "shallow-water-2d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 2
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 3

/* choices for upwinding schemes */
/*! Roe upwinding scheme */
#define _ROE_     "roe"
/*! Local Lax-Friedrich upwinding scheme */
#define _LLF_     "llf-char"

/* grid direction */
#undef _XDIR_
#undef _YDIR_
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0
#define _YDIR_ 1

/*! \def _ShallowWater2DGetFlowVar_
  Compute the flow variables (height, velocity)
  from the conserved solution vector.
  \f{equation}{
    {\bf u} = \left[\begin{array}{c} h \\ hu \\ hv \end{array}\right]
  \f}
*/
#define _ShallowWater2DGetFlowVar_(u,h,uvel,vvel) \
  { \
    h    = u[0]; \
    uvel = u[1] / h; \
    vvel = u[2] / h; \
  }

/*! \def _ShallowWater2DSetFlux_
    Set the flux vector given the flow variables
    (height, velocity).
    \f{equation}{
      {\bf F}\left({\bf u}\right)
        = \left\{\begin{array}{cc}
            \left[\begin{array}{c} hu \\ hu^2 + \frac{1}{2} gh^2 \\ huv \end{array}\right] & {\rm dir} = x \\
            \left[\begin{array}{c} hv \\ huv \\ hv^2 + \frac{1}{2} gh^2 \end{array}\right] & {\rm dir} = y \\
          \end{array}\right.
    \f}
*/
#define _ShallowWater2DSetFlux_(f,h,uvel,vvel,g,dir) \
  { \
    if (dir == _XDIR_) { \
      f[0] = (h) * (uvel); \
      f[1] = (h) * (uvel) * (uvel) + 0.5 * (g) * (h) * (h); \
      f[2] = (h) * (uvel) * (vvel); \
    } else if (dir == _YDIR_) { \
      f[0] = (h) * (vvel); \
      f[1] = (h) * (uvel) * (vvel); \
      f[2] = (h) * (vvel) * (vvel) + 0.5 * (g) * (h) * (h); \
    } \
  }

/*! \def _ShallowWater2DRoeAverage_
    Compute the Roe average of two conserved solution
    vectors.\n
*/
#define _ShallowWater2DRoeAverage_(uavg,uL,uR,p) \
  { \
    double h , uvel , vvel ; \
    double hL, uvelL, vvelL; \
    double hR, uvelR, vvelR; \
    _ShallowWater2DGetFlowVar_(uL,hL,uvelL,vvelL); \
    _ShallowWater2DGetFlowVar_(uR,hR,uvelR,vvelR); \
    h    = 0.5 * (hL    + hR   ); \
    uvel = (sqrt(hL)*uvelL + sqrt(hR)*uvelR) / (sqrt(hL) + sqrt(hR)); \
    vvel = (sqrt(hL)*vvelL + sqrt(hR)*vvelR) / (sqrt(hL) + sqrt(hR)); \
    uavg[0] = h; \
    uavg[1] = h*uvel; \
    uavg[2] = h*vvel; \
  }

/*! \def _ShallowWater2DEigenvalues_
    Compute the eigenvalues, given the conserved solution vector. The
    eigenvalues are returned as a 3x3 matrix stored in row-major format.
    It is a diagonal matrix with the eigenvalues as diagonal elements. 
    Admittedly, it is a wasteful way of storing the eigenvalues.
*/
#define _ShallowWater2DEigenvalues_(u,D,p,dir) \
  { \
    double h,uvel,vvel,c; \
    _ShallowWater2DGetFlowVar_(u,h,uvel,vvel); \
    c = sqrt(p->g*h); \
    if (dir == _XDIR_) { \
      D[0*_MODEL_NVARS_+0] = uvel-c; D[0*_MODEL_NVARS_+1] = 0;      D[0*_MODEL_NVARS_+2] = 0;      \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = uvel;   D[1*_MODEL_NVARS_+2] = 0;      \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;      D[2*_MODEL_NVARS_+2] = uvel+c; \
    } else if (dir == _YDIR_) { \
      D[0*_MODEL_NVARS_+0] = vvel-c; D[0*_MODEL_NVARS_+1] = 0;      D[0*_MODEL_NVARS_+2] = 0;      \
      D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = vvel;   D[1*_MODEL_NVARS_+2] = 0;      \
      D[2*_MODEL_NVARS_+0] = 0;      D[2*_MODEL_NVARS_+1] = 0;      D[2*_MODEL_NVARS_+2] = vvel+c; \
    } \
  }

/*! \def _ShallowWater2DLeftEigenvectors_
    Compute the matrix that has the left-eigenvectors as
    its rows. Stored in row-major format.
    Reference:
    + Ambrosi, D., "Approximation of Shallow Water Equations
      by Roe's Riemann Solver", Intl. J. Num. Meth. Fluids,
      20 (1995), pp. 157 -- 168, Page 160, Eqn. 15\n
      http://dx.doi.org/10.1002/fld.1650200205
*/
#define _ShallowWater2DLeftEigenvectors_(u,L,p,dir) \
  { \
    double h,uvel,vvel,c; \
    _ShallowWater2DGetFlowVar_(u,h,uvel,vvel); \
    c    = sqrt(p->g*h); \
    \
    if (dir == _XDIR_) { \
      L[0*_MODEL_NVARS_+0] = 0.5 + uvel/(2*c); \
      L[0*_MODEL_NVARS_+1] = -1.0/(2*c); \
      L[0*_MODEL_NVARS_+2] = 0.0; \
      \
      L[1*_MODEL_NVARS_+0] = vvel; \
      L[1*_MODEL_NVARS_+1] = 0.0; \
      L[1*_MODEL_NVARS_+2] = -1.0; \
      \
      L[2*_MODEL_NVARS_+0] = 0.5 - uvel/(2*c); \
      L[2*_MODEL_NVARS_+1] = 1.0/(2*c); \
      L[2*_MODEL_NVARS_+2] = 0.0; \
      \
    } else if (dir == _YDIR_) { \
      L[0*_MODEL_NVARS_+0] = 0.5 + vvel/(2*c); \
      L[0*_MODEL_NVARS_+1] = 0.0; \
      L[0*_MODEL_NVARS_+2] = -1.0/(2*c); \
      \
      L[1*_MODEL_NVARS_+0] = -uvel; \
      L[1*_MODEL_NVARS_+1] = 1.0; \
      L[1*_MODEL_NVARS_+2] = 0.0; \
      \
      L[2*_MODEL_NVARS_+0] = 0.5 - vvel/(2*c); \
      L[2*_MODEL_NVARS_+1] = 0.0; \
      L[2*_MODEL_NVARS_+2] = 1.0/(2*c); \
    } \
  }

/*! \def _ShallowWater2DRightEigenvectors_
    Compute the matrix that has the right-eigenvectors as
    its columns. Stored in row-major format.\n
    Reference:
    + Ambrosi, D., "Approximation of Shallow Water Equations
      by Roe's Riemann Solver", Intl. J. Num. Meth. Fluids,
      20 (1995), pp. 157 -- 168, Page 160, Eqn. 14\n
      http://dx.doi.org/10.1002/fld.1650200205
*/
#define _ShallowWater2DRightEigenvectors_(u,R,p,dir) \
  { \
    double h,uvel,vvel,c; \
    _ShallowWater2DGetFlowVar_(u,h,uvel,vvel); \
    c    = sqrt(p->g*h); \
    \
    if (dir == _XDIR_) { \
      R[0*_MODEL_NVARS_+0] = 1.0; \
      R[1*_MODEL_NVARS_+0] = uvel-c; \
      R[2*_MODEL_NVARS_+0] = vvel; \
      \
      R[0*_MODEL_NVARS_+1] = 0.0; \
      R[1*_MODEL_NVARS_+1] = 0.0; \
      R[2*_MODEL_NVARS_+1] = -1.0;\
      \
      R[0*_MODEL_NVARS_+2] = 1.0; \
      R[1*_MODEL_NVARS_+2] = uvel+c; \
      R[2*_MODEL_NVARS_+2] = vvel; \
      \
    } else if (dir == _YDIR_) { \
      R[0*_MODEL_NVARS_+0] = 1.0; \
      R[1*_MODEL_NVARS_+0] = uvel; \
      R[2*_MODEL_NVARS_+0] = vvel-c; \
      \
      R[0*_MODEL_NVARS_+1] = 0.0; \
      R[1*_MODEL_NVARS_+1] = 1.0; \
      R[2*_MODEL_NVARS_+1] = 0.0;\
      \
      R[0*_MODEL_NVARS_+2] = 1.0; \
      R[1*_MODEL_NVARS_+2] = uvel; \
      R[2*_MODEL_NVARS_+2] = vvel+c; \
    } \
  }

/*! \def ShallowWater2D
    \brief Structure containing variables and parameters specific to the 2D Shallow Water equations.
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the 2D ShallowWater equations.
*/
/*! \brief Structure containing variables and parameters specific to the 2D Shallow Water equations.
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the 2D ShallowWater equations.
*/
typedef struct shallowwater2d_parameters {
  int     bt_type;   /*!< 1 -> bottom topography is periodic, 0 -> bottom topography is not periodic */
  double  g,         /*!< Acceleration due to gravity */
          *b,        /*!< Array to store the bottom topography \f$b(x,y)\f$ */
          fhat,      /*!< Coriolis parameter */
          beta,      /*!< beta-plane approximation parameter for Coriolis force */
          D;         /*!< Channel width for Coriolis force calculation */
  char    upw_choice[_MAX_STRING_SIZE_]; /*!< Choice of upwinding scheme.\sa #_ROE_, #_LLF_*/
  /*! Function pointer to the function that computes the "upwinding" step in source term computation. To 
      understand the implementation of the gravitational source terms, see:
      + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the 
        exact conservation property for the shallow water equations", Journal of 
        Computational Physics, 208, 2005, pp. 206-227.
        http://dx.doi.org/10.1016/j.jcp.2005.02.006
  */
  int (*SourceUpwind)(double*,double*,double*,double*,int,void*,double);
} ShallowWater2D;

/*! Function to initialize the 2D ShallowWater module */
int    ShallowWater2DInitialize (void*,void*);
/*! Function to clean up the 2D ShallowWater module */
int    ShallowWater2DCleanup    (void*);

