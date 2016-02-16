/*! @file TimeExplicitRKInitialize.c
    @brief Initialize explicit Runge-Kutta time integrators
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <timeintegration.h>

/*! Initialize the explicit Runge-Kutta time integrator: Depending on the specific
    explicit Runge-Kutta (ERK) method chosen, this function allocates memory for the
    Butcher tableaux and sets their coefficients. 
*/
int TimeExplicitRKInitialize(
                              char *class,  /*!< Name of time integrator class; must match #_RK_ */
                              char *type,   /*!< Type of explicit Runge-Kutta method */
                              void *s,      /*!< Object of type #ExplicitRKParameters */
                              void *m       /*!< Unused argument */
                            )
{
  ExplicitRKParameters *params = (ExplicitRKParameters*) s;

  if (!strcmp(class,_RK_)) {
    if (!strcmp(type,_RK_1FE_)) {
      params->nstages = 1;
      params->A = (double*) calloc (params->nstages*params->nstages,sizeof(double));
      params->b = (double*) calloc (params->nstages                ,sizeof(double));
      params->c = (double*) calloc (params->nstages                ,sizeof(double));
      _ArraySetValue_(params->A,params->nstages*params->nstages,0.0);
      _ArraySetValue_(params->b,params->nstages                ,0.0);
      _ArraySetValue_(params->c,params->nstages                ,0.0);
      params->b[0] = 1.0;
    } else if (!strcmp(type,_RK_22_)) {
      params->nstages = 2;
      params->A = (double*) calloc (params->nstages*params->nstages,sizeof(double));
      params->b = (double*) calloc (params->nstages                ,sizeof(double));
      params->c = (double*) calloc (params->nstages                ,sizeof(double));
      _ArraySetValue_(params->A,params->nstages*params->nstages,0.0);
      _ArraySetValue_(params->b,params->nstages                ,0.0);
      _ArraySetValue_(params->c,params->nstages                ,0.0);
      params->A[2] = 1.0;;
      params->c[1] = 1.0;;
      params->b[0] = params->b[1] = 0.5;
    } else if (!strcmp(type,_RK_33_)) {
      params->nstages = 3;
      params->A = (double*) calloc (params->nstages*params->nstages,sizeof(double));
      params->b = (double*) calloc (params->nstages                ,sizeof(double));
      params->c = (double*) calloc (params->nstages                ,sizeof(double));
      _ArraySetValue_(params->A,params->nstages*params->nstages,0.0);
      _ArraySetValue_(params->b,params->nstages                ,0.0);
      _ArraySetValue_(params->c,params->nstages                ,0.0);
      params->A[3] = 2.0/3.0; params->A[6] = 2.0/3.0-1.0/4.0; params->A[7] = 1.0/4.0;
      params->c[1] = 2.0/3.0; params->c[2] = 2.0/3.0;
      params->b[0] = 1.0/4.0; params->b[1] = -1.0/4.0; params->b[2] = 1.0;
    } else if (!strcmp(type,_RK_44_)) {
      params->nstages = 4;
      params->A = (double*) calloc (params->nstages*params->nstages,sizeof(double));
      params->b = (double*) calloc (params->nstages                ,sizeof(double));
      params->c = (double*) calloc (params->nstages                ,sizeof(double));
      _ArraySetValue_(params->A,params->nstages*params->nstages,0.0);
      _ArraySetValue_(params->b,params->nstages                ,0.0);
      _ArraySetValue_(params->c,params->nstages                ,0.0);
      params->A[4] = 0.5; params->A[9] = 0.5; params->A[14] = 1.0;
      params->c[1] = params->c[2] = 0.5; params->c[3] = 1.0;
      params->b[0] = 1.0/6.0; params->b[1] = 1.0/3.0; params->b[2] = 1.0/3.0; params->b[3] = 1.0/6.0;
    } else if ((!strcmp(type,_RK_SSP3_)) || (!strcmp(type,_RK_TVD3_))) {
      params->nstages = 3;
      params->A = (double*) calloc (params->nstages*params->nstages,sizeof(double));
      params->b = (double*) calloc (params->nstages                ,sizeof(double));
      params->c = (double*) calloc (params->nstages                ,sizeof(double));
      _ArraySetValue_(params->A,params->nstages*params->nstages,0.0);
      _ArraySetValue_(params->b,params->nstages                ,0.0);
      _ArraySetValue_(params->c,params->nstages                ,0.0);
      params->A[3] = 1.0; params->A[6] = 0.25; params->A[7] = 0.25;
      params->c[1] = 1.0; params->c[2] = 0.5;
      params->b[0] = params->b[1] = 1.0/6.0; params->b[2] = 2.0/3.0;
    } else {
      fprintf(stderr,"Error in TimeExplicitRKInitialize(): %s is not a supported ",type);
      fprintf(stderr,"multi-stage time integration scheme of class %s.\n",class);
      return(1);
    }
  } else {
    fprintf(stderr,"Error in TimeExplicitRKInitialize(): Code should not have ");
    fprintf(stderr,"reached here for %s class of time-integrators. This is a ",class);
    fprintf(stderr,"coding mistake.\n");
    return(1);
  }
  return(0);
}
