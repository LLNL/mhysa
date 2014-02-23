/*

  Fokker-Planck Model for a 3-Bud Power System

Reference: 
+ To be added

Description:
To be added

*/

#define _FP_POWER_SYSTEM_3BUS_  "fp-power-system-3bus"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 4
#define _MODEL_NVARS_ 1

typedef struct fp_power_system__3bus_parameters {

  /* input parameters */
  double PM1, PM2,
         H1, H2,
         omegaB,
         D1, D2,
         E1, E2,
         Xd1, Xd2,
         sigma[2][2],
         lambda[2][2],
         alpha,beta;
  double *G, *B;

  /* computed/constant parameters */
  int N;
  double *Ainv;

  double pdf_integral; /* not an input          */
} FPPowerSystem3Bus;

int FPPowerSystem3BusInitialize    (void*,void*);
int FPPowerSystem3BusCleanup       (void*);
