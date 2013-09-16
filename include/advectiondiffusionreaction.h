typedef struct physics_parameters {
  double *a;                          /* advection speed along each dimension */
  int    spatial_scheme_advection;    /* choice of spatial discretization scheme for advective term */
  int    spatial_scheme_diffusion;    /* choice of spatial discretization scheme for diffusive term */
} LinearADR;

int    LinearADRInitialize        (void*,void*);
int    LinearADRCleanup           (void*,void*);
double LinearADRComputeCFL        (void*,void*,double);
double LinearADRComputeDiffNumber (void*,void*,double);
int    LinearADRAdvection         (void*,void*);
int    LinearADRDiffusion         (void*,void*);
int    LinearADRReaction          (void*,void*);

