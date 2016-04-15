/*! @file NavierStokes3DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 3D Navier Stokes equations.
*/
#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <physicalmodels/navierstokes3d.h>

/*! Compute the Roe-averaged state for the 3D Navier Stokes equations. This function 
    just calls the macro #_NavierStokes3DRoeAverage_ and is not used by any 
    functions within the 3D Navier Stokes module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes3DRoeAverage(
                              double  *uavg, /*!< The computed Roe-averaged state */
                              double  *uL,   /*!< Left state (conserved variables)*/
                              double  *uR,   /*!< Right state (conserved variables)*/
                              void    *p     /*!< Object of type #NavierStokes3D with physics-related variables */
                            )
{
  NavierStokes3D *param  = (NavierStokes3D*) p;
  _NavierStokes3DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
