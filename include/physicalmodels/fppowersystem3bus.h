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

  /* computed/constant parameters */

  double pdf_integral; /* not an input          */
} FPPowerSystem3Bus;

int FPPowerSystem3BusInitialize    (void*,void*);
int FPPowerSystem3BusCleanup       (void*);
