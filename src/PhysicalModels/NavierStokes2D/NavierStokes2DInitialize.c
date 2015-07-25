/*! @file NavierStokes2DInitialize.c
    @author Debojyoti Ghosh
    @brief Initialization of the physics-related variables and function pointers for the 2D Navier-Stokes system
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

double NavierStokes2DComputeCFL        (void*,void*,double,double);
int    NavierStokes2DFlux              (double*,double*,int,void*,double);
int    NavierStokes2DStiffFlux         (double*,double*,int,void*,double);
int    NavierStokes2DNonStiffFlux      (double*,double*,int,void*,double);
int    NavierStokes2DRoeAverage        (double*,double*,double*,void*);
int    NavierStokes2DParabolicFunction (double*,double*,void*,void*,double);
int    NavierStokes2DSource            (double*,double*,void*,void*,double);

int    NavierStokes2DJacobian          (double*,double*,void*,int,int);
int    NavierStokes2DStiffJacobian     (double*,double*,void*,int,int);

int    NavierStokes2DLeftEigenvectors  (double*,double*,void*,int);
int    NavierStokes2DRightEigenvectors (double*,double*,void*,int);

int    NavierStokes2DUpwindRoe         (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindRF          (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindLLF         (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindSWFS        (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindRusanov     (double*,double*,double*,double*,double*,double*,int,void*,double);

int    NavierStokes2DUpwinddFRoe       (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwinddFRF        (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwinddFLLF       (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwinddFRusanov   (double*,double*,double*,double*,double*,double*,int,void*,double);

int    NavierStokes2DUpwindFdFRoe      (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindFdFRF       (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindFdFLLF      (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes2DUpwindFdFRusanov  (double*,double*,double*,double*,double*,double*,int,void*,double);

int    NavierStokes2DGravityField      (void*,void*);
int    NavierStokes2DModifiedSolution  (double*,double*,int,void*,void*,double);
int    NavierStokes2DPreStep           (double*,void*,void*,double);
int    NavierStokes2DPostStage         (double*,void*,void*,double);

/*! Initialize the 2D Navier-Stokes (#NavierStokes2D) module:
    Sets the default parameters, read in and set physics-related parameters,
    and set the physics-related function pointers in #HyPar.
*/
int NavierStokes2DInitialize(
                              void *s, /*!< Solver object of type #HyPar */
                              void *m  /*!< MPI object of type #MPIVariables */
                            )
{
  HyPar           *solver  = (HyPar*)          s;
  MPIVariables    *mpi     = (MPIVariables*)   m; 
  NavierStokes2D  *physics = (NavierStokes2D*) solver->physics;
  int             ferr     = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in NavierStokes2DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in NavierStokes2DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma  = 1.4; 
  physics->Pr     = 0.72;
  physics->Re     = -1;
  physics->Minf   = 1.0;
  physics->C1     = 1.458e-6;
  physics->C2     = 110.4;
  physics->grav_x = 0.0;
  physics->grav_y = 0.0;
  physics->rho0   = 1.0;
  physics->p0     = 1.0;
  physics->HB     = 1;
  physics->R      = 1.0;
  physics->N_bv   = 0.0;
  strcpy(physics->upw_choice,"roe");

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word);                      if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word);                  if (ferr != 1) return(1);
          if (!strcmp(word, "gamma")) { 
            ferr = fscanf(in,"%lf",&physics->gamma);    if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upw_choice); if (ferr != 1) return(1);
          } else if (!strcmp(word,"Pr")) {
            ferr = fscanf(in,"%lf",&physics->Pr);       if (ferr != 1) return(1);
          } else if (!strcmp(word,"Re")) {
            ferr = fscanf(in,"%lf",&physics->Re);       if (ferr != 1) return(1);
          } else if (!strcmp(word,"Minf")) {
            ferr = fscanf(in,"%lf",&physics->Minf);     if (ferr != 1) return(1);
          } else if (!strcmp(word,"gravity")) {
            ferr = fscanf(in,"%lf",&physics->grav_x);   if (ferr != 1) return(1);
            ferr = fscanf(in,"%lf",&physics->grav_y);   if (ferr != 1) return(1);
          } else if (!strcmp(word,"rho_ref")) {
            ferr = fscanf(in,"%lf",&physics->rho0);     if (ferr != 1) return(1);
          } else if (!strcmp(word,"p_ref")) {
            ferr = fscanf(in,"%lf",&physics->p0);       if (ferr != 1) return(1);
          } else if (!strcmp(word,"HB")) {
            ferr = fscanf(in,"%d",&physics->HB);        if (ferr != 1) return(1);
            if (physics->HB==3) {
              ferr = fscanf(in,"%lf",&physics->N_bv);   if (ferr != 1) return(1);
            }
          } else if (!strcmp(word,"R")) {
            ferr = fscanf(in,"%lf",&physics->R);        if (ferr != 1) return(1);
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
  IERR MPIBroadcast_character (physics->upw_choice,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->gamma    ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Pr       ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Re       ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Minf     ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->grav_x   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->grav_y   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->rho0     ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->p0       ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->R        ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->N_bv     ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer   (&physics->HB       ,1                ,0,&mpi->world); CHECKERR(ierr);
#endif

  /* Scaling the Reynolds number with the M_inf */
  physics->Re /= physics->Minf;

  /* check that a well-balanced upwinding scheme is being used for cases with gravity */
  if (   ((physics->grav_x != 0.0) || (physics->grav_y != 0.0))
      && (strcmp(physics->upw_choice,_LLF_    )) 
      && (strcmp(physics->upw_choice,_RUSANOV_)) 
      && (strcmp(physics->upw_choice,_ROE_    ))              ) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in NavierStokes2DInitialize: %s, %s or %s upwinding is needed for flows ",_LLF_,_ROE_,_RUSANOV_);
      fprintf(stderr,"with gravitational forces.\n");
    }
    return(1);
  }
  /* check that solver has the correct choice of diffusion formulation */
  if (strcmp(solver->spatial_type_par,_NC_2STAGE_) && (physics->Re > 0)) {
    if (!mpi->rank) 
      fprintf(stderr,"Error in NavierStokes2DInitialize(): Parabolic term spatial discretization must be \"%s\"\n",_NC_2STAGE_);
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->PreStep     = NavierStokes2DPreStep;
  solver->PostStage   = NavierStokes2DPostStage;
  solver->ComputeCFL  = NavierStokes2DComputeCFL;
  solver->FFunction   = NavierStokes2DFlux;
  solver->SFunction   = NavierStokes2DSource;
  solver->UFunction   = NavierStokes2DModifiedSolution;
  if      (!strcmp(physics->upw_choice,_ROE_    )) solver->Upwind = NavierStokes2DUpwindRoe;
  else if (!strcmp(physics->upw_choice,_RF_     )) solver->Upwind = NavierStokes2DUpwindRF;
  else if (!strcmp(physics->upw_choice,_LLF_    )) solver->Upwind = NavierStokes2DUpwindLLF;
  else if (!strcmp(physics->upw_choice,_SWFS_   )) solver->Upwind = NavierStokes2DUpwindSWFS;
  else if (!strcmp(physics->upw_choice,_RUSANOV_)) solver->Upwind = NavierStokes2DUpwindRusanov;
  else {
    if (!mpi->rank) {
      fprintf(stderr,"Error in NavierStokes2DInitialize(): %s is not a valid upwinding scheme. ",
              physics->upw_choice);
      fprintf(stderr,"Choices are %s, %s, %s, %s, and %s.\n",_ROE_,_RF_,_LLF_,_SWFS_,_RUSANOV_);
    }
    return(1);
  }
  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    solver->FdFFunction = NavierStokes2DNonStiffFlux;
    solver->dFFunction  = NavierStokes2DStiffFlux;
    solver->JFunction   = NavierStokes2DStiffJacobian;
    if      (!strcmp(physics->upw_choice,_ROE_)) {
      solver->UpwinddF  = NavierStokes2DUpwinddFRoe;
      solver->UpwindFdF = NavierStokes2DUpwindFdFRoe;
    } else if (!strcmp(physics->upw_choice,_RF_)) {
      solver->UpwinddF  = NavierStokes2DUpwinddFRF;
      solver->UpwindFdF = NavierStokes2DUpwindFdFRF;
    } else if (!strcmp(physics->upw_choice,_LLF_)) {
      solver->UpwinddF  = NavierStokes2DUpwinddFLLF;
      solver->UpwindFdF = NavierStokes2DUpwindFdFLLF;
    } else if (!strcmp(physics->upw_choice,_RUSANOV_)) {
      solver->UpwinddF  = NavierStokes2DUpwinddFRusanov;
      solver->UpwindFdF = NavierStokes2DUpwindFdFRusanov;
    } else {
      if (!mpi->rank) {
        fprintf(stderr,"Error in NavierStokes2DInitialize(): %s is not a valid upwinding scheme ",
                physics->upw_choice);
        fprintf(stderr,"for use with split hyperbolic flux form. Use %s, %s, %s, or %s.\n",
                _ROE_,_RF_,_LLF_,_RUSANOV_);
      }
      return(1);
    }
  } else solver->JFunction      = NavierStokes2DJacobian;
  solver->AveragingFunction     = NavierStokes2DRoeAverage;
  solver->GetLeftEigenvectors   = NavierStokes2DLeftEigenvectors;
  solver->GetRightEigenvectors  = NavierStokes2DRightEigenvectors;

  /* set the value of gamma in all the boundary objects */
  int n;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  for (n = 0; n < solver->nBoundaryZones; n++)  boundary[n].gamma = physics->gamma;

  /* hijack the main solver's dissipation function pointer
   * to this model's own function, since it's difficult to express 
   * the dissipation terms in the general form                      */
  solver->ParabolicFunction = NavierStokes2DParabolicFunction;

  /* allocate array to hold the gravity field */
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int d, size = 1; for (d=0; d<_MODEL_NDIMS_; d++) size *= (dim[d] + 2*ghosts);
  physics->grav_field_f = (double*) calloc (size, sizeof(double));
  physics->grav_field_g = (double*) calloc (size, sizeof(double));
  physics->fast_jac     = (double*) calloc (2*size*_MODEL_NVARS_*_MODEL_NVARS_,sizeof(double));
  physics->solution     = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));
  IERR NavierStokes2DGravityField(solver,mpi); CHECKERR(ierr);

  return(0);
}
