#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

double FPDoubleWellComputeCFL        (void*,void*,double);
double FPDoubleWellComputeDiffNumber (void*,void*,double);
int    FPDoubleWellAdvection         (double*,double*,int,void*,double);
int    FPDoubleWellDiffusion         (double*,double*,int,void*,double);
int    FPDoubleWellUpwind            (double*,double*,double*,double*,int,void*);

int FPDoubleWellInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  FPDoubleWell  *physics = (FPDoubleWell*)  solver->physics;
  int           ierr     = 0,i;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in FPDoubleWellInitializeO(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in FPDoubleWellInitializeO(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* reading physical model specific inputs - all processes */
  FILE *in;
  if (!mpi->rank) printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) {
    fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
    char word[_MAX_STRING_SIZE_];
    ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
        if (!strcmp(word, "q")) {
          /* read diffusion coefficient */
          ierr = fscanf(in,"%lf",&physics->q);
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
  solver->ComputeCFL         = FPDoubleWellComputeCFL;
  solver->ComputeDiffNumber  = FPDoubleWellComputeDiffNumber;
  solver->FFunction          = FPDoubleWellAdvection;
  solver->GFunction          = FPDoubleWellDiffusion;
  solver->SFunction          = NULL;
  solver->Upwind             = FPDoubleWellUpwind;

  return(0);
}
