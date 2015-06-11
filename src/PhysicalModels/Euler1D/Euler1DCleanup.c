/*! @file Euler1DCleanup.c
    @author Debojyoti Ghosh
    @brief Contains the function to clean up the 1D Euler physics module
*/

#include <stdlib.h>
#include <physicalmodels/euler1d.h>

/*! Function to clean up all physics-related allocations for the 1D Euler equations */
int Euler1DCleanup(
                   void *s /*!< Solver object of type #HyPar */
                  )
{
  Euler1D *param  = (Euler1D*) s;
  free(param->grav_field);
  free(param->fast_jac);
  free(param->solution);
  return(0);
}
