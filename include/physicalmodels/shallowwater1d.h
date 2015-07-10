/*! @file shallowwater1d.h
    @brief 1D Shallow Water Equations
    @author Debojyoti Ghosh

    1D Shallow Water Equations\n

    \f{equation}{
      \frac {\partial} {\partial t} \left[\begin{array}{c} h \\ hu \end{array} \right]
    + \frac {\partial} {\partial x} \left[\begin{array}{c} hu \\ hu^2 + \frac{1}{2}gh^2  \end{array} \right] 
    = \left[\begin{array}{c} 0 \\ -ghb_x \end{array}\right]
    \f}
    where \f$h\f$ is the height, \f$u\f$ is the velocity, \f$b(x)\f$ is the bottom 
    topography, and \f$g\f$ is the gravitational constant.\n\n

    For the treatment of the source term (well-balanced formulation), refer to:\n
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the 
      exact conservation property for the shallow water equations", Journal of 
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/

#include <basic.h>
#include <math.h>
#include <matops.h>

/*! \def _SHALLOW_WATER_1D_
    1D Shallow Water equations
*/
#define _SHALLOW_WATER_1D_ "shallow-water-1d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 1
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 2

/* choices for upwinding schemes */
/*! Roe upwinding scheme */
#define _ROE_     "roe"
/*! Local Lax-Friedrich upwinding scheme */
#define _LLF_     "llf-char"

/* grid direction */
#undef _XDIR_
/*! dimension corresponding to the \a x spatial dimension */
#define _XDIR_ 0

/*! \def _ShallowWater1DGetFlowVar_
  Compute the flow variables (height, velocity)
  from the conserved solution vector.
*/
#define _ShallowWater1DGetFlowVar_(u,h,v) \
  { \
    h = u[0]; \
    v = u[1] / h; \
  }

/*! \def _ShallowWater1DSetFlux_
    Set the flux vector given the flow variables
    (height, velocity).
*/
#define _ShallowWater1DSetFlux_(f,h,v,g) \
  { \
    f[0] = (h) * (v); \
    f[1] = (h) * (v) * (v) + 0.5 * (g) * (h) * (h); \
  }

/*! \def _ShallowWater1DRoeAverage_
    Compute the Roe average of two conserved solution
    vectors.\n
    Note: Implemented as the arithmetic average right now.
*/
#define _ShallowWater1DRoeAverage_(uavg,uL,uR,p) \
  { \
    double h ,v; \
    double hL,vL; \
    double hR,vR; \
    hL = uL[0]; \
    vL = uL[1] / hL; \
    hR = uR[0]; \
    vR = uR[1] / hR; \
    h  = 0.5 * (hL + hR); \
    v  = 0.5 * (vL + vR); \
    uavg[0] = h; \
    uavg[1] = h*v; \
  }

/*! \def _ShallowWater1DEigenvalues_
    Compute the eigenvalues, given the conserved solution vector. The
    eigenvalues are returned as a 3x3 matrix stored in row-major format.
    It is a diagonal matrix with the eigenvalues as diagonal elements. 
    Admittedly, it is a wasteful way of storing the eigenvalues.
*/
#define _ShallowWater1DEigenvalues_(u,D,p,dir) \
  { \
    double h,v,c; \
    h = u[0]; \
    v = u[1] / h; \
    c = sqrt(p->g*h); \
    D[0*_MODEL_NVARS_+0] = (v+c);  D[0*_MODEL_NVARS_+1] = 0;     \
    D[1*_MODEL_NVARS_+0] = 0;      D[1*_MODEL_NVARS_+1] = (v-c); \
  }

/*! \def _ShallowWater1DLeftEigenvectors_
    Compute the matrix that has the left-eigenvectors as
    its rows. Stored in row-major format.
*/
#define _ShallowWater1DLeftEigenvectors_(u,L,p,dir) \
  { \
    double R[_MODEL_NVARS_*_MODEL_NVARS_]; \
    _ShallowWater1DRightEigenvectors_(u,R,p,dir); \
    _MatrixInvert_(R,L,_MODEL_NVARS_); \
  }

/*! \def _ShallowWater1DRightEigenvectors_
    Compute the matrix that has the right-eigenvectors as
    its columns. Stored in row-major format.
*/
#define _ShallowWater1DRightEigenvectors_(u,R,p,dir) \
  { \
    double h,v,c; \
    h = u[0]; \
    v = u[1] / h; \
    c    = sqrt(p->g*h); \
    R[0*_MODEL_NVARS_+0] = 1.0; \
    R[1*_MODEL_NVARS_+0] = v+c; \
    R[0*_MODEL_NVARS_+1] = 1.0; \
    R[1*_MODEL_NVARS_+1] = v-c; \
  }

/*! \def ShallowWater1D
    \brief Structure containing variables and parameters specific to the 1D Shallow Water equations.
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the 1D ShallowWater equations.
*/
/*! \brief Structure containing variables and parameters specific to the 1D Shallow Water equations.
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the 1D ShallowWater equations.
*/
typedef struct shallowwater1d_parameters {
  double  g;         /*!< Acceleration due to gravity */
  int     bt_type;   /*!< 1 -> bottom topography is periodic, 0 -> bottom topography is not periodic */
  double  *b;        /*!< Array to store the bottom topography \f$b(x)\f$ */
  char    upw_choice[_MAX_STRING_SIZE_]; /*!< Choice of upwinding scheme.\sa #_ROE_, #_LLF_*/
  /*! Function pointer to the function that computes the "upwinding" step in source term computation. To 
      understand the implementation of the gravitational source terms, see:
      + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the 
        exact conservation property for the shallow water equations", Journal of 
        Computational Physics, 208, 2005, pp. 206-227.
        http://dx.doi.org/10.1016/j.jcp.2005.02.006
  */
  int (*SourceUpwind)(double*,double*,double*,double*,int,void*,double);
} ShallowWater1D;

/*! Function to initialize the 1D ShallowWater module */
int    ShallowWater1DInitialize (void*,void*);
/*! Function to clean up the 1D ShallowWater module */
int    ShallowWater1DCleanup    (void*);

