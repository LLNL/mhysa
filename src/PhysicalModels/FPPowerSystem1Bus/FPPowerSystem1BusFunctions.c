#include <math.h>
#include <physicalmodels/fppowersystem1bus.h>

double FPPowerSystem1BusDriftFunction(int dir,void *p,double x,double y, double t)
{
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) p;

  double drift = 0;
  if (dir == _XDIR_) {
    drift = params->omegaB * (y - params->omegaS);
  } else if (dir == _YDIR_) {
    drift = (params->omegaS/(2*params->H)) 
          * (params->Pm_avg - params->Pmax*sin(x) - params->D*(y-params->omegaS));
  }

  return drift;
}

double FPPowerSystem1BusDissipationFunction(int dir1,int dir2,void *p,double t)
{
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) p;

  double dissp = 0;
  if (dir1 == _YDIR_) {
    if (dir2 == _XDIR_) {
      dissp = (params->sigma*params->sigma*params->omegaS*params->omegaS)/(4*params->H*params->H) 
              * params->lambda * params->lambda * params->omegaB;
    } else if (dir2 == _YDIR_) {
      dissp = (params->sigma*params->sigma*params->omegaS*params->omegaS)/(4*params->H*params->H) 
              * params->lambda * (1.0 - params->lambda * (params->D*params->omegaS)/(2*params->H));
    }
  }

  return(dissp);
}
