#include <math.h>
#include <basic.h>
#include <physicalmodels/navierstokes2d.h>

int NavierStokes2DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  NavierStokes2D *param  = (NavierStokes2D*) p;
  _NavierStokes2DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
