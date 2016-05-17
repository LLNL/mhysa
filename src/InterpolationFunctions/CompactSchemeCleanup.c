/*! @file CompactSchemeCleanup.c
    @brief Cleans up allocations specific to compact finite-difference methods
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Cleans up all allocations related to the compact finite difference methods.
*/
int CompactSchemeCleanup(void *s /*!< CompactScheme object of type #CompactScheme */ )
{
  CompactScheme *compact   = (CompactScheme*) s;

  if (compact->A) free(compact->A);
  if (compact->B) free(compact->B);
  if (compact->C) free(compact->C);
  if (compact->R) free(compact->R);

  if (compact->sendbuf) free(compact->sendbuf);
  if (compact->recvbuf) free(compact->recvbuf);

  return(0);
}
