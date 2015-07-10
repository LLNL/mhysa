#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

double FPDoubleWellComputeCFL        (void*,void*,double,double);
double FPDoubleWellComputeDiffNumber (void*,void*,double,double);
int    FPDoubleWellAdvection         (double*,double*,int,void*,double);
int    FPDoubleWellDiffusion         (double*,double*,int,void*,double);
int    FPDoubleWellUpwind            (double*,double*,double*,double*,
                                      double*,double*,int,void*,double);
int    FPDoubleWellPostStep          (double*,void*,void*,double,int);
int    FPDoubleWellPrintStep         (void*,void*,double);

int FPDoubleWellInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  FPDoubleWell  *physics = (FPDoubleWell*)  solver->physics;
  int           ferr     = 0;
  _DECLARE_IERR_;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in FPDoubleWellInitializeO(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in FPDoubleWellInitializeO(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) {
      fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
      return(1);
    } else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "q")) {
            /* read diffusion coefficient */
            ferr = fscanf(in,"%lf",&physics->q);
            if (ferr != 1) return(1);
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
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
  }

#ifndef serial
  IERR MPIBroadcast_double(&physics->q,1,0,&mpi->world); CHECKERR(ierr);
#endif

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in FPDoubleWellInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = FPDoubleWellComputeCFL;
  solver->ComputeDiffNumber  = FPDoubleWellComputeDiffNumber;
  solver->FFunction          = FPDoubleWellAdvection;
  solver->GFunction          = FPDoubleWellDiffusion;
  solver->Upwind             = FPDoubleWellUpwind;
  solver->PostStep           = FPDoubleWellPostStep;
  solver->PrintStep          = FPDoubleWellPrintStep;

  /* Calculate and print the PDF integral of the initial solution */
  IERR FPDoubleWellPostStep(solver->u,solver,mpi,0.0,0);CHECKERR(ierr);
  IERR FPDoubleWellPrintStep(solver,mpi,0.0);           CHECKERR(ierr);
  
  return(0);
}
