#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int MUSCLInitialize(void *s,void *m)
{
  MUSCLParameters  *muscl   = (MUSCLParameters*) s;
  /* hard coding these parameters for now */
  /* modify to read from an input file later */
  muscl->eps = 1e-3;
  return(0);
}
