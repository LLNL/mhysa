/*! @file BCCleanup.c
    @author Debojyoti Ghosh
    @brief Function to clean up boundary-conditions related variables
*/

#include <stdlib.h>
#include <boundaryconditions.h>

/*! Cleans up a boundary object of type #DomainBoundary */
int BCCleanup(void *b /*!< Boundary object of type #DomainBoundary*/)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  free(boundary->xmin);
  free(boundary->xmax);
  free(boundary->is);
  free(boundary->ie);
  if (boundary->DirichletValue) free(boundary->DirichletValue);
  if (boundary->SpongeValue   ) free(boundary->SpongeValue   );
  if (boundary->FlowDensity  ) free(boundary->FlowDensity  );
  if (boundary->FlowVelocity  ) free(boundary->FlowVelocity  );
  if (boundary->UnsteadyDirichletSize) free(boundary->UnsteadyDirichletSize);
  if (boundary->UnsteadyDirichletData) free(boundary->UnsteadyDirichletData);
  return(0);
}
