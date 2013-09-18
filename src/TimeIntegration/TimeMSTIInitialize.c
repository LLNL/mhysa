#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <timeintegration.h>

int TimeMSTIInitialize(char *class,char *type,void *s)
{
  MSTIParameters *params = (MSTIParameters*) s;

  if (!strcmp(class,_RK_)) {
    if (!strcmp(type,_RK_1FE_)) {
      params->nstages = 1;
      params->A = (double*) calloc (params->nstages*params->nstages,sizeof(double));
      params->b = (double*) calloc (params->nstages                ,sizeof(double));
      params->A[0] = 0.0;
      params->b[0] = 1.0;
    } else {
      fprintf(stderr,"Error in TimeMSTIInitialize(): %s is not a supported ",type);
      fprintf(stderr,"multi-stage time integration scheme of class %s.\n",class);
      return(1);
    }
  } else {
    fprintf(stderr,"Error in TimeMSTIInitialize(): %s is not a supported ",class);
    fprintf(stderr,"multi-stage time integrator class.\n");
    return(1);
  }
  return(0);
}
