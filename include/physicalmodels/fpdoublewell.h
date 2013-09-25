/*

Reference: 
+ Miller, R.N., Carter E.F., Blue S.T., "Data assimilation into
  nonlinear stochastic models", Tellus (1999), 51 A, 167-194

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
