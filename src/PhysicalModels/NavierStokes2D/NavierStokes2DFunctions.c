/*! @file NavierStokes2DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 2D Navier Stokes equations.
*/
#include <math.h>
#include <basic.h>
#include <physicalmodels/navierstokes2d.h>

/*! Compute the Roe-averaged state for the 2D Navier Stokes equations. This function 
    just calls the macro #_NavierStokes2DRoeAverage_ and is not used by any 
    functions within the 2D Navier Stokes module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes2DRoeAverage(
                              double  *uavg, /*!< The computed Roe-averaged state */
                              double  *uL,   /*!< Left state (conserved variables)*/
                              double  *uR,   /*!< Right state (conserved variables)*/
                              void    *p     /*!< Object of type #NavierStokes2D with physics-related variables */
                            )
{
  NavierStokes2D *param  = (NavierStokes2D*) p;
  _NavierStokes2DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
