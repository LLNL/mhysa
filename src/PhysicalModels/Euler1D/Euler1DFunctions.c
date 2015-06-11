/*! @file Euler1DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 1D Euler equations.
*/

#include <math.h>
#include <basic.h>
#include <physicalmodels/euler1d.h>

/*! Compute the Roe-averaged state for the 1D Euler equations. This function 
    just calls the macro #_Euler1DRoeAverage_ and is not used by any 
    functions within the 1D Euler module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int Euler1DRoeAverage(
                      double  *uavg, /*!< The computed Roe-averaged state */
                      double  *uL,   /*!< Left state (conserved variables)*/
                      double  *uR,   /*!< Right state (conserved variables)*/
                      void    *p     /*!< Object of type #Euler1D with physics-related variables */
                     )
{
  Euler1D *param  = (Euler1D*) p;
  _Euler1DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
