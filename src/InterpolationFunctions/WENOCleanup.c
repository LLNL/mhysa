/*! @file WENOCleanup.c
    @brief Cleans up allocations specific to WENO-type methods
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Cleans up all allocations related to the WENO-type methods.
*/
int WENOCleanup(void *s /*!< WENO object of type #WENOParameters */ )
{
  WENOParameters  *weno   = (WENOParameters*) s;

  if (weno->A) free(weno->A);
  if (weno->B) free(weno->B);
  if (weno->C) free(weno->C);
  if (weno->R) free(weno->R);

  if (weno->sendbuf) free(weno->sendbuf);
  if (weno->recvbuf) free(weno->recvbuf);

  free(weno->offset);
  free(weno->w1);
  free(weno->w2);
  free(weno->w3);

  return(0);
}
