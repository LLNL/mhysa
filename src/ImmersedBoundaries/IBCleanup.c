/*! @file IBCleanup.c
    @brief Clean up immersed boundaries-related allocations.
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <immersedboundaries.h>

/*!
*/
int IBCleanup(void *s /*!< Object of type #ImmersedBoundary */)
{
  ImmersedBoundary *ib = (ImmersedBoundary*) s;
  if (!ib) return(0);

  free(ib->body->surface);
  free(ib->body);

  if (ib->n_boundary_nodes > 0) free(ib->boundary);
  if (ib->nfacets_local > 0) free(ib->fmap);

  return(0);
}
