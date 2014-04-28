#include <stdlib.h>
#include <boundaryconditions.h>

int BCCleanup(void *b)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  free(boundary->xmin);
  free(boundary->xmax);
  free(boundary->is);
  free(boundary->ie);
  if (boundary->DirichletValue) free(boundary->DirichletValue);
  if (boundary->SpongeValue   ) free(boundary->SpongeValue   );
  if (boundary->FlowVelocity  ) free(boundary->FlowVelocity  );
  if (boundary->UnsteadyDirichletSize) free(boundary->UnsteadyDirichletSize);
  if (boundary->UnsteadyDirichletData) free(boundary->UnsteadyDirichletData);
  return(0);
}
