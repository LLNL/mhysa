#include <stdio.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

int FPPowerSystem3BusPrintStep(void* s,void *m,double t)
{
  HyPar               *solver = (HyPar*)              s;
  FPPowerSystem3Bus   *params = (FPPowerSystem3Bus*)  solver->physics;
  printf("Domain integral of the probability density function: %1.16E\n",params->pdf_integral);
  return(0);
}
