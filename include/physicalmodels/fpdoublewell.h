/*
  Fokker-Planck Model for Double Well

Reference: 
+ Miller, R.N., Carter E.F., Blue S.T., "Data assimilation into
  nonlinear stochastic models", Tellus (1999), 51 A, 167-194

  dP     d[f(x)P]     1   d^2 P
  --  =  --------  +  - q -----
  dt        dx        2   dx^2

  f(x) = 4x(x^2-1) (drift)
  q = constant (input)


*/

#define _FP_DOUBLE_WELL_  "fp-double-well"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 1
#define _MODEL_NVARS_ 1

typedef struct fp_double_well_parameters {
  double q;            /* diffusion coefficient */
  double pdf_integral; /* not an input          */
} FPDoubleWell;

#define drift(x) (4.0*(x)*(1.0-(x)*(x))) 

int FPDoubleWellInitialize    (void*,void*);
int FPDoubleWellCleanup       (void*);
