#define _LINEAR_ADVECTION_DIFFUSION_REACTION_  "linear-advection-diffusion-reaction"

typedef struct linearadr_parameters {
  double *a;  /* advection speed for each variable along each dimension        */
  double *d;  /* diffusion coefficient for each variable along  each dimension */

} LinearADR;

int    LinearADRInitialize        (void*,void*);
int    LinearADRCleanup           (void*);

