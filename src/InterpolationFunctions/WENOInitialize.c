#include <stdlib.h>
#include <string.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int WENOInitialize(void *s,void *m, char *scheme)
{
  HyPar           *solver = (HyPar*) s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;

  /* hard coding these parameters for now */
  /* modify to read from an input file later */
  weno->mapped      = 1;
  weno->borges      = 0;
  weno->yc          = 0;
  weno->no_limiting = 0;
  weno->eps         = 1e-6;
  weno->p           = 2.0;

  if (!strcmp(scheme,_FIFTH_ORDER_CRWENO_)) {
    int size = 1, d;
    for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+1);
    size *= solver->nvars;

    weno->A = (double*) calloc (size, sizeof(double));
    weno->B = (double*) calloc (size, sizeof(double));
    weno->C = (double*) calloc (size, sizeof(double));
    weno->R = (double*) calloc (size, sizeof(double));

    weno->sendbuf = (double*) calloc (size, sizeof(double));
    weno->recvbuf = (double*) calloc (size, sizeof(double));

  } else {

    weno->A = weno->B = weno->C = weno->R = NULL;
    weno->sendbuf = weno->recvbuf = NULL;

  }

  return(0);
}
