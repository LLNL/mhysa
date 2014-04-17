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
  if (boundary->FlowVelocity  ) free(boundary->FlowVelocity  );
  return(0);
}
