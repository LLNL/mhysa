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

  double term = (params->sigma*params->sigma*params->omegaS*params->omegaS) / (4.0*params->H*params->H);
  double expterm = exp(-t/params->lambda);
  double dissp = 0;
  if (dir1 == _YDIR_) {
    if (dir2 == _XDIR_) {
      /* dissp = 0.0; */
      dissp = term * (params->lambda*params->omegaB) * (params->lambda*(1-expterm) - t*expterm);
      /* dissp = term * params->lambda*params->omegaB*params->lambda; */
    } else if (dir2 == _YDIR_) {
      double gamma = params->D * params->omegaS / (2.0*params->H);
      /* dissp =   term * (params->lambda/(params->lambda*gamma+1.0))
              * (1.0 - exp(-(gamma+1.0/params->lambda)*t)); */
      dissp = term * (params->lambda*(1-expterm) 
              + ( gamma*params->lambda * (t*expterm-params->lambda*(1-expterm)) ) );
      /* dissp = term * params->lambda * (1 - params->lambda*gamma); */
    }
  }

  return(dissp);
}
