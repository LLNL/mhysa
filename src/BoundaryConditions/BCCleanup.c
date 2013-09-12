#include <stdlib.h>
#include <boundaryconditions.h>

int BCCleanup(void *b)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  free(boundary->xmin);
  free(boundary->xmax);
  free(boundary->is);
  free(boundary->ie);

  return(0);
}
