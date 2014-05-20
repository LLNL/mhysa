#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <physicalmodels/euler2d.h>
#include <mpivars.h>
#include <hypar.h>

double Euler2DComputeCFL        (void*,void*,double,double);
int    Euler2DFlux              (double*,double*,int,void*,double);
int    Euler2DUpwindRoe         (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler2DUpwindRF          (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler2DUpwindLLF         (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler2DUpwindSWFS        (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler2DRoeAverage        (double*,double*,double*,void*);
int    Euler2DLeftEigenvectors  (double*,double*,void*,int);
int    Euler2DRightEigenvectors (double*,double*,void*,int);

int Euler2DInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)          s;
  MPIVariables  *mpi     = (MPIVariables*)   m; 
  Euler2D       *physics = (Euler2D*) solver->physics;
  int           ferr     = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in Euler2DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in Euler2DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma = 1.4; 
  strcpy(physics->upw_choice,"roe");

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "gamma")) { 
            ferr = fscanf(in,"%lf",&physics->gamma); if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upw_choice); if (ferr != 1) return(1);
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
  IERR MPIBroadcast_double    (&physics->gamma,1,0,&mpi->world);                      CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->upw_choice,_MAX_STRING_SIZE_,0,&mpi->world);  CHECKERR(ierr);
#endif

  /* initializing physical model-specific functions */
  solver->ComputeCFL  = Euler2DComputeCFL;
  solver->FFunction   = Euler2DFlux;
  if      (!strcmp(physics->upw_choice,_ROE_ )) solver->Upwind = Euler2DUpwindRoe;
  else if (!strcmp(physics->upw_choice,_RF_  )) solver->Upwind = Euler2DUpwindRF;
  else if (!strcmp(physics->upw_choice,_LLF_ )) solver->Upwind = Euler2DUpwindLLF;
  else if (!strcmp(physics->upw_choice,_SWFS_)) solver->Upwind = Euler2DUpwindSWFS;
  else {
    fprintf(stderr,"Error in Euler2DInitialize(): %s is not a valid upwinding scheme.\n",
            physics->upw_choice);
    return(1);
  }
  solver->AveragingFunction     = Euler2DRoeAverage;
  solver->GetLeftEigenvectors   = Euler2DLeftEigenvectors;
  solver->GetRightEigenvectors  = Euler2DRightEigenvectors;

  /* set the value of gamma in all the boundary objects */
  int n;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  for (n = 0; n < solver->nBoundaryZones; n++)  boundary[n].gamma = physics->gamma;

  return(0);
}
