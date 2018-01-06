/*! @file NavierStokes3DCleanup.c
    @author Debojyoti Ghosh
    @brief Clean up the 3D Navier Stokes module
*/
#include <stdlib.h>
#include <physicalmodels/navierstokes3d.h>

/*! Function to clean up all allocations in the 3D Navier
    Stokes module.
*/
int NavierStokes3DCleanup(void *s /*!< Object of type #NavierStokes3D*/)
{
  NavierStokes3D  *param  = (NavierStokes3D*) s;

  return(0);
}
