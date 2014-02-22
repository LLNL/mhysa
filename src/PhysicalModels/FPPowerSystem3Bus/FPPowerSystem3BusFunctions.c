#include <math.h>
#include <physicalmodels/fppowersystem3bus.h>

int FPPowerSystem3BusDriftFunction(int dir,void *p,double *x, double t, double *drift)
{
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;

  int N = params->N;
  double theta1, Omega1, F1, Tau1, gamma1, S1, VR1, VI1;
  double theta2, Omega2, F2, Tau2, gamma2, S2, VR2, VI2;

  theta1 = x[0];
  theta2 = x[1];
  Omega1 = x[2];
  Omega2 = x[3];

  VR1 = params->Ainv[0*N+0]*params->E1*sin(theta1) - params->Ainv[0*N+1]*params->E1*cos(theta1) 
      + params->Ainv[0*N+2]*params->E2*sin(theta2) - params->Ainv[0*N+3]*params->E2*cos(theta2);
  VI1 = params->Ainv[1*N+0]*params->E1*sin(theta1) - params->Ainv[1*N+1]*params->E1*cos(theta1) 
      + params->Ainv[1*N+2]*params->E2*sin(theta2) - params->Ainv[1*N+3]*params->E2*cos(theta2);
  VR2 = params->Ainv[2*N+0]*params->E1*sin(theta1) - params->Ainv[2*N+1]*params->E1*cos(theta1) 
      + params->Ainv[2*N+2]*params->E2*sin(theta2) - params->Ainv[2*N+3]*params->E2*cos(theta2);
  VI2 = params->Ainv[3*N+0]*params->E1*sin(theta1) - params->Ainv[3*N+1]*params->E1*cos(theta1) 
      + params->Ainv[3*N+2]*params->E2*sin(theta2) - params->Ainv[3*N+3]*params->E2*cos(theta2);

  F1 = params->PM1 / (2.0*params->H1*params->omegaB);
  F2 = params->PM2 / (2.0*params->H2*params->omegaB);

  Tau1 = 0.0;
  Tau2 = 0.0;

  gamma1 = params->D1 / (2*params->H1*params->omegaB);
  gamma2 = params->D2 / (2*params->H2*params->omegaB);

  S1 = params->E1/(2.0*params->H1*params->omegaB*params->Xd1) * (VR1*sin(theta1) - VI1*cos(theta1));
  S2 = params->E2/(2.0*params->H2*params->omegaB*params->Xd2) * (VR2*sin(theta2) - VI2*cos(theta2));

  drift[0] = Omega1;
  drift[1] = Omega2;
  drift[2] = F1 + Tau1 - gamma1*Omega1 - S1;
  drift[3] = F2 + Tau2 - gamma2*Omega2 - S2;

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
