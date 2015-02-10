#include <stdio.h>
#include <string.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeGetAuxSolutions(int *N, double **uaux, void *s,int n)
{
  HyPar           *solver = (HyPar*) s;
  TimeIntegration *TS     = (TimeIntegration*) solver->time_integrator;

  if (n >= 0) {
    if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
      GLMGEEParameters *params = (GLMGEEParameters*) solver->msti;
      *uaux = TS->U[params->r+n];
    }
  } else {
    if (!TS) *N = 0;
    else {
      if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
        GLMGEEParameters *params = (GLMGEEParameters*) solver->msti;
        *N = params->r - 1;
      } else *N = 0;
    }
  }

  return(0);
}
