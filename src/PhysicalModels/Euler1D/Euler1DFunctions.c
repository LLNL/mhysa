#include <math.h>
#include <basic.h>
#include <physicalmodels/euler1d.h>

int Euler1DRoeAverage(double *uavg,double *uL,double *uR,void *p)
{
  Euler1D *param  = (Euler1D*) p;
  _Euler1DRoeAverage_(uavg,uL,uR,param);
  return(0);
}
