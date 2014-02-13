#include <math.h>
#include <physicalmodels/fppowersystem3bus.h>

int FPPowerSystem3BusDriftFunction(int dir,void *p,double *x, double t, double *drift)
{
//  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;

  drift[0] = 0;
  drift[1] = 0;
  drift[2] = 0;
  drift[3] = 0;

  return(0);
}

int FPPowerSystem3BusDissipationFunction(int dir,void *p,double t, double *dissp)
{
//  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;

  dissp[0] = 0;
  dissp[1] = 0;
  dissp[2] = 0;
  dissp[3] = 0;

  return(0);
}
