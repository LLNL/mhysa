#include <timeintegration.h>

int TimeStep(void *ts)
{
  TimeIntegration *TS  = (TimeIntegration*) ts;
  int             ierr = 0;
//  if (TS->TimeIntegrate) ierr = TS->TimeIntegrate(TS);
  return(ierr);
}

