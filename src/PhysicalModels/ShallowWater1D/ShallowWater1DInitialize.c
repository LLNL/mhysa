/*! @file ShallowWater1DInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize the 1D shallow water equations module.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <mpivars.h>
#include <hypar.h>

double ShallowWater1DComputeCFL (void*,void*,double,double);
int    ShallowWater1DFlux       (double*,double*,int,void*,double);
int    ShallowWater1DSource     (double*,double*,void*,void*,double);
int    ShallowWater1DUpwindRoe  (double*,double*,double*,double*,double*,double*,int,void*,double);
int    ShallowWater1DUpwindLLF  (double*,double*,double*,double*,double*,double*,int,void*,double);

int    ShallowWater1DJacobian          (double*,double*,void*,int,int);
int    ShallowWater1DRoeAverage        (double*,double*,double*,void*);
int    ShallowWater1DLeftEigenvectors  (double*,double*,void*,int);
int    ShallowWater1DRightEigenvectors (double*,double*,void*,int);

int    ShallowWater1DTopography        (void*,void*);
int    ShallowWater1DSourceUpwindLLF   (double*,double*,double*,double*,int,void*,double);
int    ShallowWater1DSourceUpwindRoe   (double*,double*,double*,double*,int,void*,double);

int    ShallowWater1DModifiedSolution  (double*,double*,int,void*,void*,double);
int    ShallowWater1DWriteTopography   (void*,void*);

/*! Function to initialize the 1D shallow water equations (#ShallowWater1D) module: 
    Sets the default parameters, read in and set physics-related parameters, 
    and set the physics-related function pointers in #HyPar.
*/
int ShallowWater1DInitialize(
                      void *s, /*!< Solver object of type #HyPar */
                      void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                     )
{
  HyPar           *solver  = (HyPar*)         s;
  MPIVariables    *mpi     = (MPIVariables*)  m; 
  ShallowWater1D  *physics = (ShallowWater1D*)       solver->physics;
  int             ferr, d;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in ShallowWater1DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in ShallowWater1DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->g       = 1.0; 
  physics->bt_type = 0;
  strcpy(physics->upw_choice,"roe");

  /* reading physical model specific inputs */
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
          if (!strcmp(word, "gravity")) { 
            ferr = fscanf(in,"%lf",&physics->g); 
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "topography_type")) { 
            ferr = fscanf(in,"%d",&physics->bt_type); 
            if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upw_choice); 
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
  IERR MPIBroadcast_double    (&physics->g        ,1,0,&mpi->world);                  CHECKERR(ierr);
  IERR MPIBroadcast_integer   (&physics->bt_type  ,1,0,&mpi->world);                  CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->upw_choice,_MAX_STRING_SIZE_,0,&mpi->world);  CHECKERR(ierr);
#endif

  /* initializing physical model-specific functions */
  solver->ComputeCFL = ShallowWater1DComputeCFL;
  solver->FFunction  = ShallowWater1DFlux;
  solver->SFunction  = ShallowWater1DSource;
  solver->UFunction  = ShallowWater1DModifiedSolution;
  solver->JFunction  = ShallowWater1DJacobian;
  if      (!strcmp(physics->upw_choice,_ROE_ )) solver->Upwind = ShallowWater1DUpwindRoe;
  else if (!strcmp(physics->upw_choice,_LLF_ )) solver->Upwind = ShallowWater1DUpwindLLF;
  else {
    if (!mpi->rank) fprintf(stderr,"Error in ShallowWater1DInitialize(): %s is not a valid upwinding scheme.\n",
                            physics->upw_choice);
    return(1);
  }
  solver->AveragingFunction     = ShallowWater1DRoeAverage;
  solver->GetLeftEigenvectors   = ShallowWater1DLeftEigenvectors;
  solver->GetRightEigenvectors  = ShallowWater1DRightEigenvectors;
  solver->PhysicsOutput         = ShallowWater1DWriteTopography;
   
  if      (!strcmp(physics->upw_choice,_LLF_ )) physics->SourceUpwind = ShallowWater1DSourceUpwindLLF;
  else if (!strcmp(physics->upw_choice,_ROE_ )) physics->SourceUpwind = ShallowWater1DSourceUpwindRoe;

  /* allocate array to hold the bottom topography field */
  physics->b = (double*) calloc (solver->npoints_local_wghosts, sizeof(double));
  IERR ShallowWater1DTopography(solver,mpi); CHECKERR(ierr);

  return(0);
}
