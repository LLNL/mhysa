/*! @file ShallowWater1DCleanup.c
    @author Debojyoti Ghosh
    @brief Contains the function to clean up the 1D shallow water physics module
*/

#include <stdlib.h>
#include <physicalmodels/shallowwater1d.h>

/*! Function to clean up all physics-related allocations for the 1D shallow water equations */
int ShallowWater1DCleanup(
                   void *s /*!< Solver object of type #HyPar */
                  )
{
  ShallowWater1D *param  = (ShallowWater1D*) s;
  free(param->b);
  return(0);
}
