/*! @file ShallowWater2DFlux.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the hyperbolic flux for the 2D shallow water equations over the domain.
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <hypar.h>

/*! Compute the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}\left({\bf u}\right) 
      = \left\{\begin{array}{cc} 
          \left[\begin{array}{c} hu \\ hu^2 + \frac{1}{2} gh^2 \\ huv \end{array}\right] & {\rm dir} = x \\
          \left[\begin{array}{c} hv \\ huv \\ hv^2 + \frac{1}{2} gh^2 \end{array}\right] & {\rm dir} = y \\
        \end{array}\right.
    \f}
*/
int ShallowWater2DFlux(
                        double  *f, /*!< Array to hold the computed flux (same size and layout as u) */
                        double  *u, /*!< Array containing the conserved solution */
                        int     dir,/*!< Spatial dimension */
                        void    *s, /*!< Solver object of type #HyPar */
                        double  t   /*!< Current time */
                      )
{
  HyPar             *solver = (HyPar*)   s;
  ShallowWater2D    *param  = (ShallowWater2D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  nvars   = _MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double h, uvel, vvel;
    _ShallowWater2DGetFlowVar_((u+nvars*p),h,uvel,vvel);
    _ShallowWater2DSetFlux_((f+nvars*p),h,uvel,vvel,param->g,dir);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
