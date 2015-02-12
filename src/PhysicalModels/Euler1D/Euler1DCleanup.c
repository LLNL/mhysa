#include <stdlib.h>
#include <physicalmodels/euler1d.h>

int Euler1DCleanup(void *s)
{
  Euler1D *param  = (Euler1D*) s;
  free(param->grav_field);
  free(param->fast_jac);
  free(param->solution);
  return(0);
}
