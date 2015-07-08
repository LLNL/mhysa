/*! @file ShallowWater1DFlux.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the hyperbolic flux for the 1D shallow water equations over the domain.
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <hypar.h>

/*! Compute the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}\left({\bf u}\right) = \left[\begin{array}{c} hu \\ hu^2 + \frac{1}{2} gh^2 \end{array}\right]
    \f}
*/
int ShallowWater1DFlux(
                        double  *f, /*!< Array to hold the computed flux (same size and layout as u) */
                        double  *u, /*!< Array containing the conserved solution */
                        int     dir,/*!< Spatial dimension (unused since this is a 1D system) */
                        void    *s, /*!< Solver object of type #HyPar */
                        double  t   /*!< Current time */
                      )
{
  HyPar             *solver = (HyPar*)   s;
  ShallowWater1D    *param  = (ShallowWater1D*) solver->physics;
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
    double h, v;
    _ShallowWater1DGetFlowVar_((u+nvars*p),h,v);
    _ShallowWater1DSetFlux_((f+nvars*p),h,v,param->g);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
