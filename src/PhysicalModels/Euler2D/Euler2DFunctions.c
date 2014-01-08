#include <math.h>
#include <basic.h>
#include <physicalmodels/euler2d.h>

int Euler2DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  Euler2D *param  = (Euler2D*) p;
  _Euler2DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
