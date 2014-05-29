#include <string.h>
#include <basic.h>
#include <interpolation.h>
#include <hypar.h>

/*
  For non-linear schemes (like WENO,CRWENO,HCWENO), precalculate
  the non-linear interpolation coefficients (ie weights)
*/

int InterpSetLimiterVar(double *fC,double *u,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  _DECLARE_IERR_;

  if (  (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_WENO_  ))
      ||(!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_))
      ||(!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_)) ) {
    IERR weno->CalculateWENOWeights(fC,u,x,dir,s,m);
  }
  return(0);
}
