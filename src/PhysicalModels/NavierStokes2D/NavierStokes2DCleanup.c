/*! @file NavierStokes2DCleanup.c
    @author Debojyoti Ghosh
    @brief Clean up the 2D Navier Stokes module
*/
#include <stdlib.h>
#include <physicalmodels/navierstokes2d.h>

/*! Function to clean up all allocations in the 2D Navier
    Stokes module.
*/
int NavierStokes2DCleanup(void *s /*!< Object of type #NavierStokes2D*/)
{
  NavierStokes2D  *param  = (NavierStokes2D*) s;

  free(param->grav_field_f);
  free(param->grav_field_g);
  free(param->fast_jac);
  free(param->solution);
  return(0);
}
