#include <math.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d.h>

int Numa2DCalculateStandardAtmosphere_1(void *p,double z,double *ExnerP,double *P,double *rho,double *T)
{
  Numa2D *physics = (Numa2D*) p;

  double R      = physics->R;
  double gamma  = physics->gamma;
  double g      = physics->g;

  /* reference quantities at zero altitude */
  double P0, T0; 
  P0   = physics->Pref;
  T0   = physics->Tref;

  double inv_gamma_m1 = 1.0/(gamma-1.0);
  double Cp = gamma * inv_gamma_m1 * R;

  double theta  = T0;
  *ExnerP = 1.0 - (g/(Cp*theta))*z;
  *P      = P0 * raiseto((*ExnerP),gamma*inv_gamma_m1);
  *rho    = (P0/(R*theta)) * raiseto((*ExnerP),inv_gamma_m1);
  *T      = (*rho) * theta;

  return(0);
}

int Numa2DCalculateStandardAtmosphere_2(void *p,double z,double *ExnerP,double *P,double *rho,double *T)
{
  Numa2D *physics = (Numa2D*) p;

  double R      = physics->R;
  double gamma  = physics->gamma;
  double g      = physics->g;

  /* reference quantities at zero altitude */
  double P0, T0; 
  P0 = physics->Pref;
  T0 = physics->Tref;

  double BV = 0.01; /* Brunt-Vaisala frequency */
  double inv_gamma_m1 = 1.0/(gamma-1.0);
  double Cp = gamma * inv_gamma_m1 * R;

  double term   = BV*BV*z/g;
  double theta  = T0 * exp(term);
  *ExnerP = 1.0 + (g*g/(Cp*T0*BV*BV)) * (exp(-term) - 1.0);
  *P      = P0 * raiseto((*ExnerP),gamma*inv_gamma_m1);
  *rho    = (P0/(R*theta)) * raiseto((*ExnerP),inv_gamma_m1);
  *T      = (*rho) * theta;

  return(0);
}
