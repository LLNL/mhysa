/*

  1D Euler Equations for Inviscid, Compressible Flows

    
  d   [ rho   ]   d   [   rho*u    ]
  --  [ rho*u ] + --  [rho*u*u + p ] = 0
  dt  [   e   ]   dx  [ (e+p)*u    ]

  Equation of state:
           p         1
    e = -------  +   - rho * u^2
        gamma-1      2

  Choices for upwinding:
  "roe"       Roe upwinding
  "rf-char"   Roe-fixed
  "llf-char"  Local Lax-Friedrich

*/

#include <basic.h>

#define _EULER_1D_  "euler1d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 1
#define _MODEL_NVARS_ 3

typedef struct euler1d_parameters {
  double  gamma;  /* Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
} Euler1D;

int    Euler1DInitialize (void*,void*);
int    Euler1DCleanup    (void*);

