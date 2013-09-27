/*
  Linear Advection-Diffusion-Reaction

Reference:
+ Hundsdorfer & Verwer, "Numerical Solution of Time Dependent Advection-
  Diffusion-Reaction Euqations", Springer-Verlag Berlin Heidelberg, 2010

  du         du          d^2 u
  --  + a_i ----  = nu_i ------   + k u
  dt        dx_i         dx_i^2

*/


#define _LINEAR_ADVECTION_DIFFUSION_REACTION_  "linear-advection-diffusion-reaction"

typedef struct linearadr_parameters {
  double *a;  /* advection speed for each variable along each dimension        */
  double *d;  /* diffusion coefficient for each variable along  each dimension */

} LinearADR;

int    LinearADRInitialize        (void*,void*);
int    LinearADRCleanup           (void*);

