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

  VR1 = params->Ainv[0*N+0]*params->E1*sin(theta1)/params->Xd1 - params->Ainv[0*N+1]*params->E1*cos(theta1)/params->Xd1 
      + params->Ainv[0*N+2]*params->E2*sin(theta2)/params->Xd2 - params->Ainv[0*N+3]*params->E2*cos(theta2)/params->Xd2;
  VI1 = params->Ainv[1*N+0]*params->E1*sin(theta1)/params->Xd1 - params->Ainv[1*N+1]*params->E1*cos(theta1)/params->Xd1 
      + params->Ainv[1*N+2]*params->E2*sin(theta2)/params->Xd2 - params->Ainv[1*N+3]*params->E2*cos(theta2)/params->Xd2;
  VR2 = params->Ainv[2*N+0]*params->E1*sin(theta1)/params->Xd1 - params->Ainv[2*N+1]*params->E1*cos(theta1)/params->Xd1 
      + params->Ainv[2*N+2]*params->E2*sin(theta2)/params->Xd2 - params->Ainv[2*N+3]*params->E2*cos(theta2)/params->Xd2;
  VI2 = params->Ainv[3*N+0]*params->E1*sin(theta1)/params->Xd1 - params->Ainv[3*N+1]*params->E1*cos(theta1)/params->Xd1 
      + params->Ainv[3*N+2]*params->E2*sin(theta2)/params->Xd2 - params->Ainv[3*N+3]*params->E2*cos(theta2)/params->Xd2;

  F1 = params->PM1 / (2.0*params->H1);
  F2 = params->PM2 / (2.0*params->H2);

  Tau1 = 0.0;
  Tau2 = 0.0;

  gamma1 = params->D1 / (2*params->H1);
  gamma2 = params->D2 / (2*params->H2);

  S1 = params->E1/(2.0*params->H1*params->Xd1) * (VR1*sin(theta1) - VI1*cos(theta1));
  S2 = params->E2/(2.0*params->H2*params->Xd2) * (VR2*sin(theta2) - VI2*cos(theta2));

  drift[0] = params->omegaB * (Omega1-1.0);
  drift[1] = params->omegaB * (Omega2-1.0);
  drift[2] = F1 + Tau1 - gamma1*(Omega1-1.0) - S1;
  drift[3] = F2 + Tau2 - gamma2*(Omega2-1.0) - S2;

  return(0);
}

int FPPowerSystem3BusDissipationFunction(int dir1,int dir2,void *p,double t, double *dissp)
{
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;
  int i; for (i=0; i<_MODEL_NDIMS_*_MODEL_NDIMS_; i++) dissp[i] = 0.0;

  double sigma11 = params->sigma[0][0];
  double sigma12 = params->sigma[0][1];
  double sigma21 = params->sigma[1][0];
  double sigma22 = params->sigma[1][1];

  double lambda11 = params->lambda[0][0];
  double lambda12 = params->lambda[0][1];
  double lambda21 = params->lambda[1][0];
  double lambda22 = params->lambda[1][1];

  double gamma1 = params->D1 / (2*params->H1*params->omegaB);
  double gamma2 = params->D2 / (2*params->H2*params->omegaB);

  double term11 = 1.0 - exp((lambda11*t)/(1-lambda11*(gamma1+gamma2)));
  double term12 = 1.0 - exp((lambda12*t)/(1-lambda12*(gamma1+gamma2)));
  double term21 = 1.0 - exp((lambda21*t)/(1-lambda21*(gamma1+gamma2)));
  double term22 = 1.0 - exp((lambda22*t)/(1-lambda22*(gamma1+gamma2)));

  dissp[2*_MODEL_NDIMS_+2] = (sigma11*sigma11*lambda11)/(1-lambda11*(gamma1+gamma2))*term11;
  dissp[2*_MODEL_NDIMS_+3] = (sigma12*sigma12*lambda12)/(1-lambda12*(gamma1+gamma2))*term12;
  dissp[3*_MODEL_NDIMS_+2] = (sigma21*sigma21*lambda21)/(1-lambda21*(gamma1+gamma2))*term21;
  dissp[3*_MODEL_NDIMS_+3] = (sigma22*sigma22*lambda22)/(1-lambda22*(gamma1+gamma2))*term22;

  return(0);
}
