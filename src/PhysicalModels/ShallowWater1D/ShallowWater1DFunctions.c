/*! @file ShallowWater1DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 1D shallow water equations.
*/

#include <math.h>
#include <basic.h>
#include <physicalmodels/shallowwater1d.h>

/*! Compute the Roe-averaged state for the 1D shallow water equations. This function 
    just calls the macro #_ShallowWater1DRoeAverage_ and is not used by any 
    functions within the 1D shallow water module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int ShallowWater1DRoeAverage(
                      double  *uavg, /*!< The computed Roe-averaged state */
                      double  *uL,   /*!< Left state (conserved variables)*/
                      double  *uR,   /*!< Right state (conserved variables)*/
                      void    *p     /*!< Object of type #ShallowWater1D with physics-related variables */
                     )
{
  ShallowWater1D *param  = (ShallowWater1D*) p;
  _ShallowWater1DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
