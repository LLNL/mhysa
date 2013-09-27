#include <stdio.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

int FPPowerSystemPrintStep(void* s,void *m,double t)
{
  HyPar          *solver = (HyPar*)         s;
  FPPowerSystem  *params = (FPPowerSystem*) solver->physics;
  printf("Domain integral of the probability density function: %1.16E\n",params->pdf_integral);
  return(0);
}
