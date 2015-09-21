/*! @file FPPowerSystem3BusCleanup.c
    @author Debojyoti Ghosh
    @brief Function to clean up allocations for the #FPPowerSystem3Bus system
*/

#include <stdlib.h>
#include <physicalmodels/fppowersystem3bus.h>

/*! Clean up allocations in the #FPPowerSystem3Bus system */
int FPPowerSystem3BusCleanup(
                              void *s /*!< Object of type #FPPowerSystem3Bus */
                            )
{
  FPPowerSystem3Bus *physics = (FPPowerSystem3Bus*) s;
  free(physics->G);
  free(physics->B);
  return(0);
}
