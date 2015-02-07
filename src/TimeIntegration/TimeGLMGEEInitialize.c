#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <timeintegration.h>

int TimeGLMGEEInitialize(char *class,char *type,void *s,void *m)
{
  GLMGEEParameters *params = (GLMGEEParameters*) s;
  MPIVariables     *mpi    = (MPIVariables*) m;
  int i,j;

  if (!strcmp(class,_GLM_GEE_)) {

    if (!strcmp(type,_GLM_GEE_23_)) {
      params->nstages = 3;
      params->r       = 2;
    } else {
      fprintf(stderr,"Error in TimeGLMGEEInitialize(): %s is not a supported ",type);
      fprintf(stderr,"multi-stage time integration scheme of class %s.\n",class);
      return(1);
    }
  
    int s = params->nstages;
    int r = params->r;

    params->A_yyt = (double*) calloc (s*s,sizeof(double));
    params->B_yyt = (double*) calloc (s*r,sizeof(double));
    params->C_yyt = (double*) calloc (s*r,sizeof(double));
    params->D_yyt = (double*) calloc (r*r,sizeof(double));
    params->c_yyt = (double*) calloc (s  ,sizeof(double));
    _ArraySetValue_(params->A_yyt,s*s,0.0);
    _ArraySetValue_(params->B_yyt,s*r,0.0);
    _ArraySetValue_(params->C_yyt,s*r,0.0);
    _ArraySetValue_(params->D_yyt,r*r,0.0);
    _ArraySetValue_(params->c_yyt,s  ,0.0);
    
    params->A_yeps = (double*) calloc (s*s,sizeof(double));
    params->B_yeps = (double*) calloc (s*r,sizeof(double));
    params->C_yeps = (double*) calloc (s*r,sizeof(double));
    params->D_yeps = (double*) calloc (r*r,sizeof(double));
    params->c_yeps = (double*) calloc (s  ,sizeof(double));
    _ArraySetValue_(params->A_yeps,s*s,0.0);
    _ArraySetValue_(params->B_yeps,s*r,0.0);
    _ArraySetValue_(params->C_yeps,s*r,0.0);
    _ArraySetValue_(params->D_yeps,r*r,0.0);
    _ArraySetValue_(params->c_yeps,s  ,0.0);
    
    if (!strcmp(type,_GLM_GEE_23_)) {

      params->A_yeps[1*s+0] = 1.0;
      params->A_yeps[2*s+0] = params->A_yeps[2*s+1] = 0.25;

      params->B_yeps[0*s+0] = params->B_yeps[0*s+1] = 1.0/12.0; params->B_yeps[0*s+2] =  5.0/6.0;
      params->B_yeps[1*s+0] = params->B_yeps[1*s+1] = 1.0/12.0; params->B_yeps[1*s+2] = -1.0/6.0;

      params->C_yeps[0*r+0] = 1.0;
      params->C_yeps[1*r+0] = 1.0;
      params->C_yeps[2*r+0] = 1.0;  params->C_yeps[2*r+1] = -1.0;

      params->D_yeps[0*r+0] = 1.0;
      params->D_yeps[1*r+1] = 1.0;

      params->A_yyt[1*s+0] = 1.0;
      params->A_yyt[2*s+0] = params->A_yyt[2*s+1] = 0.25;

      params->B_yyt[0*s+0] = params->B_yyt[0*s+1] = 1.0/12.0; params->B_yyt[0*s+2] =  5.0/6.0;
      params->B_yyt[1*s+0] = params->B_yyt[1*s+1] = 1.0/6.0;  params->B_yyt[1*s+2] =  2.0/3.0;

      params->C_yyt[0*r+0] =  1.0;
      params->C_yyt[1*r+0] = -9.0;  params->C_yyt[1*r+1] = 10.0;
      params->C_yyt[2*r+0] =  2.0;  params->C_yyt[2*r+1] = -1.0;

      params->D_yyt[0*r+0] = 1.0;
      params->D_yyt[1*r+1] = 1.0;

      params->gamma = 0.0;
    }

    for (i=0; i<s; i++) {
      for (j=0; j<s; j++) {
        params->c_yyt[i]  += params->A_yyt [i*s+j];
        params->c_yeps[i] += params->A_yeps[i*s+j];
      }
    }

    if (!mpi->rank) {
      FILE *in;
      int ferr;
      in = fopen("glm_gee.inp","r");
      strcat(params->ee_mode,_GLM_GEE_YEPS_);
      if (in) {
        printf("Reading GLM-GEE method parameters from glm_gee.inp.\n");
        char word[_MAX_STRING_SIZE_];
        ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
        if (!strcmp(word,"begin")) {
          while (strcmp(word,"end")) {
            ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
            if (!strcmp(word,"ee_mode")) { 
              ferr = fscanf(in,"%s",params->ee_mode); 
              if (ferr != 1) return(1); 
            } else if (strcmp(word,"end")) {
              char useless[_MAX_STRING_SIZE_];
              ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
              printf("Warning: keyword %s in file \"glm_gee.inp\" with value %s not ",word,useless);
              printf("recognized or extraneous. Ignoring.\n");
            }
          }
        } else {
          fprintf(stderr,"Error: Illegal format in file \"glm_gee.inp\".\n");
          return(1);
        }
        fclose(in);
        if (strcmp(params->ee_mode,_GLM_GEE_YEPS_) && strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
          fprintf(stderr,"Error in TimeGLMGEEInitialize(): %s is not a valid value for ",params->ee_mode);
          fprintf(stderr,"ee_mode. Acceptable inputs are %s or %s.\n",_GLM_GEE_YEPS_,_GLM_GEE_YYT_);
          strcat(params->ee_mode,_GLM_GEE_YEPS_);
        }
      }
    }
    IERR MPIBroadcast_character(params->ee_mode,_MAX_STRING_SIZE_,0,&mpi->world);

    if (!strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
      params->A = params->A_yyt;
      params->B = params->B_yyt;
      params->C = params->C_yyt;
      params->D = params->D_yyt;
      params->c = params->c_yyt;
    } else {
      params->A = params->A_yeps;
      params->B = params->B_yeps;
      params->C = params->C_yeps;
      params->D = params->D_yeps;
      params->c = params->c_yeps;
    }

  } else {
    fprintf(stderr,"Error in TimeGLMGEEInitialize(): Code should not have ");
    fprintf(stderr,"reached here for %s class of time-integrators. This is a ",class);
    fprintf(stderr,"coding mistake.\n");
    return(1);
  }
  return(0);
}
