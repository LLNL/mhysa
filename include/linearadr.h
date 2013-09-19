typedef struct physics_parameters {
  double *a;                          /* advection speed along each dimension */

} LinearADR;

int    LinearADRInitialize        (void*,void*);
int    LinearADRCleanup           (void*);

