/*! @file TimeStep.c
    @brief Advance one time step
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <timeintegration.h>

/*!
  Advance one time step.
*/
int TimeStep(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration *TS  = (TimeIntegration*) ts;
  _DECLARE_IERR_;
  if (TS->TimeIntegrate) { IERR TS->TimeIntegrate(TS); CHECKERR(ierr); }
  return(0);
}
