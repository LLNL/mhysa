#include <string.h>
#include <interpolation.h>

/*
  For non-linear schemes (like WENO,CRWENO,HCWENO), set
  the variable on whose basis the limiting will be carried out.
  NULL means use the variable being interpolated (the usual).
*/

int InterpSetLimiterVar(void *s,char *scheme,double *var)
{
  if (  (!strcmp(scheme,_FIFTH_ORDER_WENO_  ))
      ||(!strcmp(scheme,_FIFTH_ORDER_CRWENO_))
      ||(!strcmp(scheme,_FIFTH_ORDER_HCWENO_)) ) {
    WENOParameters  *weno = (WENOParameters*) s;
    weno->var = var;
  }
  return(0);
}
