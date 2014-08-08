#include <basic.h>
#include <timeintegration.h>

int TimeStep(void *ts)
{
  TimeIntegration *TS  = (TimeIntegration*) ts;
  _DECLARE_IERR_;
  if (TS->TimeIntegrate) { IERR TS->TimeIntegrate(TS); CHECKERR(ierr); }
  return(0);
}
