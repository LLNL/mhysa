#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int WENOInitialize(void *s,void *m)
{
  WENOParameters  *weno   = (WENOParameters*) s;

  /* hard coding these parameters for now */
  /* modify to read from an input file later */
  weno->mapped      = 0;
  weno->borges      = 0;
  weno->yc          = 0;
  weno->no_limiting = 0;
  weno->eps         = 1e-6;
  weno->p           = 2.0;

  return(0);
}
