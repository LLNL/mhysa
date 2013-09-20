#define _FP_DOUBLE_WELL_  "fp-double-well"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 1
#define _MODEL_NVARS_ 1

typedef struct fp_double_well_parameters {
  double q;
} FPDoubleWell;

int FPDoubleWellInitialize    (void*,void*);
int FPDoubleWellCleanup       (void*);
