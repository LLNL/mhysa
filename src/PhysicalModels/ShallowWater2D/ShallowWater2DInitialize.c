/*! @file ShallowWater2DInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize the 2D shallow water equations module.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <mpivars.h>
#include <hypar.h>

double ShallowWater2DComputeCFL (void*,void*,double,double);
int    ShallowWater2DFlux       (double*,double*,int,void*,double);
int    ShallowWater2DSource     (double*,double*,void*,void*,double);
int    ShallowWater2DUpwindRoe  (double*,double*,double*,double*,double*,double*,int,void*,double);
int    ShallowWater2DUpwindLLF  (double*,double*,double*,double*,double*,double*,int,void*,double);

int    ShallowWater2DJacobian          (double*,double*,void*,int,int);
int    ShallowWater2DRoeAverage        (double*,double*,double*,void*);
int    ShallowWater2DLeftEigenvectors  (double*,double*,void*,int);
int    ShallowWater2DRightEigenvectors (double*,double*,void*,int);

int    ShallowWater2DTopography        (void*,void*);
int    ShallowWater2DSourceUpwindLLF   (double*,double*,double*,double*,int,void*,double);
int    ShallowWater2DSourceUpwindRoe   (double*,double*,double*,double*,int,void*,double);

int    ShallowWater2DModifiedSolution  (double*,double*,int,void*,void*,double);
int    ShallowWater2DWriteTopography   (void*,void*);

/*! Function to initialize the 2D shallow water equations (#ShallowWater2D) module: 
    Sets the default parameters, read in and set physics-related parameters, 
    and set the physics-related function pointers in #HyPar.
*/
int ShallowWater2DInitialize(
                      void *s, /*!< Solver object of type #HyPar */
                      void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                     )
{
  HyPar           *solver  = (HyPar*)         s;
  MPIVariables    *mpi     = (MPIVariables*)  m; 
  ShallowWater2D  *physics = (ShallowWater2D*)       solver->physics;
  int             ferr, d;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in ShallowWater2DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in ShallowWater2DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
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
  solver->ComputeCFL = ShallowWater2DComputeCFL;
  solver->FFunction  = ShallowWater2DFlux;
  solver->SFunction  = ShallowWater2DSource;
  solver->UFunction  = ShallowWater2DModifiedSolution;
  solver->JFunction  = ShallowWater2DJacobian;
  if      (!strcmp(physics->upw_choice,_ROE_ )) solver->Upwind = ShallowWater2DUpwindRoe;
  else if (!strcmp(physics->upw_choice,_LLF_ )) solver->Upwind = ShallowWater2DUpwindLLF;
  else {
    if (!mpi->rank) fprintf(stderr,"Error in ShallowWater2DInitialize(): %s is not a valid upwinding scheme.\n",
                            physics->upw_choice);
    return(1);
  }
  solver->AveragingFunction     = ShallowWater2DRoeAverage;
  solver->GetLeftEigenvectors   = ShallowWater2DLeftEigenvectors;
  solver->GetRightEigenvectors  = ShallowWater2DRightEigenvectors;
  solver->PhysicsOutput         = ShallowWater2DWriteTopography;
   
  if      (!strcmp(physics->upw_choice,_LLF_ )) physics->SourceUpwind = ShallowWater2DSourceUpwindLLF;
  else if (!strcmp(physics->upw_choice,_ROE_ )) physics->SourceUpwind = ShallowWater2DSourceUpwindRoe;

  /* allocate array to hold the bottom topography field */
  physics->b = (double*) calloc (solver->npoints_local_wghosts, sizeof(double));
  IERR ShallowWater2DTopography(solver,mpi); CHECKERR(ierr);

  return(0);
}
