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

  double sigma  = params->sigma;
  double omegaS = params->omegaS;
  double omegaB = params->omegaB;
  double lambda = params->lambda;
  double H      = params->H;
  double D      = params->D;

  double dissp = 0;
  if (dir1 == _YDIR_) {

    double term = (sigma*sigma*omegaS*omegaS) / (4.0*H*H);
    double expterm = exp(-t/lambda);

    if (dir2 == _XDIR_) {
      dissp = term * (lambda*omegaB) * (lambda*(1-expterm) - t*expterm);
      /* dissp = term * lambda*omegaB*lambda; */
    } else if (dir2 == _YDIR_) {
      double gamma = D*omegaS / (2.0*H);
      dissp = term * (lambda*(1-expterm) + (gamma*lambda*(t*expterm-lambda*(1-expterm))) );
      /* dissp = term * params->lambda * (1 - params->lambda*gamma); */
    }

  }

  return(dissp);
}
