#include <stdio.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

int FPPowerSystem1BusPrintStep(void* s,void *m,double t)
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*)  solver->physics;
  printf("Domain integral of the probability density function: %1.16E\n",params->pdf_integral);
  return(0);
}
