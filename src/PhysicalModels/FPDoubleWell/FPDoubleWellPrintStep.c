#include <stdio.h>
#include <physicalmodels/fpdoublewell.h>
#include <hypar.h>

int FPDoubleWellPrintStep(void* s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  FPDoubleWell  *params = (FPDoubleWell*) solver->physics;
  printf("Domain integral of the probability density function: %1.16E\n",params->pdf_integral);
  return(0);
}
