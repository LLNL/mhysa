#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int WENOCleanup(void *s)
{
  WENOParameters  *weno   = (WENOParameters*) s;

  /* hard coding these parameters for now */
  /* modify to read from an input file later */
  weno->mapped      = 1;
  weno->borges      = 0;
  weno->yc          = 0;
  weno->no_limiting = 0;
  weno->eps         = 1e-6;
  weno->p           = 2.0;

  if (weno->A) free(weno->A);
  if (weno->B) free(weno->B);
  if (weno->C) free(weno->C);
  if (weno->R) free(weno->R);

  if (weno->sendbuf) free(weno->sendbuf);
  if (weno->recvbuf) free(weno->recvbuf);

  return(0);
}
