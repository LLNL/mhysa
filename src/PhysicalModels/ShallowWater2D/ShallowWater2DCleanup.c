/*! @file ShallowWater2DCleanup.c
    @author Debojyoti Ghosh
    @brief Contains the function to clean up the 2D shallow water physics module
*/

#include <stdlib.h>
#include <physicalmodels/shallowwater2d.h>

/*! Function to clean up all physics-related allocations for the 2D shallow water equations */
int ShallowWater2DCleanup(
                   void *s /*!< Solver object of type #HyPar */
                  )
{
  ShallowWater2D *param  = (ShallowWater2D*) s;
  free(param->b);
  return(0);
}
