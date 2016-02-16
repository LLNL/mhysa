/*! @file TimeExplicitRKCleanup.c
    @brief Clean up explicit Runge-Kutta time integrator
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <string.h>
#include <timeintegration.h>

/*! Clean up allocations related to explicit Runge-Kutta time integration:
    This function frees up the arrays for the Butcher tableaux. */
int TimeExplicitRKCleanup(void *s /*!< Object of type #ExplicitRKParameters*/ )
{
  ExplicitRKParameters *params = (ExplicitRKParameters*) s;
  if (params->A) free(params->A);
  if (params->b) free(params->b);
  if (params->c) free(params->c);
  return(0);
}
