#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

double Euler1DComputeCFL (void*,void*,double,double);
int    Euler1DFlux       (double*,double*,int,void*,double);
int    Euler1DUpwindRoe  (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler1DUpwindRF   (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler1DUpwindLLF  (double*,double*,double*,double*,double*,double*,int,void*,double);

int Euler1DInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  Euler1D       *physics = (Euler1D*)       solver->physics;
  int           ierr     = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in Euler1DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in Euler1DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma = 1.4; 
  strcpy(physics->upw_choice,"roe");


  /* reading physical model specific inputs - all processes */
  FILE *in;
  if (!mpi->rank) printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) {
    printf("Warning: File \"physics.inp\" not found. Using default values.\n");
  } else {
    char word[_MAX_STRING_SIZE_];
    ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
        if (!strcmp(word, "gamma")) { 
          ierr = fscanf(in,"%lf",&physics->gamma); 
          if (ierr != 1) return(1);
        } else if (!strcmp(word,"upwinding")) {
          ierr = fscanf(in,"%s",physics->upw_choice); 
          if (ierr != 1) return(1);
        } else if (strcmp(word,"end")) {
          char useless[_MAX_STRING_SIZE_];
          ierr = fscanf(in,"%s",useless); if (ierr != 1) return(ierr);
          printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",word,useless);
          printf("recognized or extraneous. Ignoring.\n");
        }
      }
	  } else {
    	fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
  }
  fclose(in);

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = Euler1DComputeCFL;
  solver->FFunction          = Euler1DFlux;
  if      (!strcmp(physics->upw_choice,"roe"))      solver->Upwind = Euler1DUpwindRoe;
  else if (!strcmp(physics->upw_choice,"rf-char"))  solver->Upwind = Euler1DUpwindRF;
  else if (!strcmp(physics->upw_choice,"llf-char")) solver->Upwind = Euler1DUpwindLLF;
  else {
    fprintf(stderr,"Error in Euler1DInitialize(): %s is not a valid upwinding scheme.\n",
            physics->upw_choice);
    return(1);
  }

  return(0);
}
