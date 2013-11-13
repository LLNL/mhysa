#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <physicalmodels/navierstokes3d.h>

int NavierStokes3DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  NavierStokes3D *param  = (NavierStokes3D*) p;
  _NavierStokes3DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
