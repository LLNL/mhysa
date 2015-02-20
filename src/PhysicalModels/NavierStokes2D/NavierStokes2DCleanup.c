#include <stdlib.h>
#include <physicalmodels/navierstokes2d.h>

int NavierStokes2DCleanup(void *s)
{
  NavierStokes2D  *param  = (NavierStokes2D*) s;

  free(param->grav_field_f);
  free(param->grav_field_g);
  free(param->fast_jac);
  free(param->solution);
  return(0);
}
