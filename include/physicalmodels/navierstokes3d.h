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

#define _NAVIER_STOKES_3D_  "navierstokes3d"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 3
#define _MODEL_NVARS_ 5

/* choices for upwinding schemes */
#define _ROE_   "roe"
#define _RF_    "rf-char"
#define _LLF_   "llf-char"

/* grid directions */
#define _XDIR_ 0
#define _YDIR_ 1
#define _ZDIR_ 2

typedef struct navierstokes3d_parameters {
  double  gamma;  /* Ratio of heat capacities */
  char    upw_choice[_MAX_STRING_SIZE_]; /* choice of upwinding */
} NavierStokes3D;

int    NavierStokes3DInitialize (void*,void*);
int    NavierStokes3DCleanup    (void*);

