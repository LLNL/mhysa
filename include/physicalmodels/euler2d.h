/*

  1D Euler Equations for Inviscid, Compressible Flows

    
  d   [ rho   ]   d   [   rho*u    ]    d  [   rho*v   ]
  --  [ rho*u ] + --  [rho*u*u + p ] +  -- [  rho*u*v  ] = 0
  dt  [ rho*v ]   dx  [   rho*u*v  ]    dy [rho*v*v + p]
      [   e   ]       [   (e+p)*u  ]       [  (e+p)v   ]

  Equation of state:
           p         1
    e = -------  +   - rho * (u^2 + v^2)
        gamma-1      2

  Choices for upwinding:
  "roe"       Roe upwinding
  "rf-char"   Roe-fixed
  "llf-char"  Local Lax-Friedrich

*/

#include <basic.h>

#define _EULER_2D_  "euler2d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 2
#define _MODEL_NVARS_ 4

/* choices for upwinding schemes */
#define _ROE_   "roe"
#define _RF_    "rf-char"
#define _LLF_   "llf-char"

/* directions */
#define _XDIR_ 0
#define _YDIR_ 1

typedef struct euler2d_parameters {
  double  gamma;  /* Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
} Euler2D;

int    Euler2DInitialize (void*,void*);
int    Euler2DCleanup    (void*);

