typedef struct physics_parameters {
  double *a;  /* advection speed for each variable along each dimension        */
  double *d;  /* diffusion coefficient for each variable along  each dimension */

} LinearADR;

int    LinearADRInitialize        (void*,void*);
int    LinearADRCleanup           (void*);

