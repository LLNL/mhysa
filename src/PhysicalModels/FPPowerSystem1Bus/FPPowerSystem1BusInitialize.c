#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystem1BusComputeCFL        (void*,void*,double,double);
double FPPowerSystem1BusComputeDiffNumber (void*,void*,double,double);
int    FPPowerSystem1BusAdvection         (double*,double*,int,void*,double);
int    FPPowerSystem1BusDiffusion         (double*,double*,int,int,void*,double);
int    FPPowerSystem1BusUpwind            (double*,double*,double*,double*,
                                           double*,double*,int,void*,double);

int FPPowerSystem1BusInitialize(void *s,void *m)
{
  HyPar               *solver  = (HyPar*)             s;
  MPIVariables        *mpi     = (MPIVariables*)      m; 
  FPPowerSystem1Bus   *physics = (FPPowerSystem1Bus*) solver->physics;
  int                 ferr;
  _DECLARE_IERR_;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in FPPowerSystem1BusInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in FPPowerSystem1BusInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values of model parameters */
  physics->omegaS = 1.0;
  physics->omegaB = 120 * (4.0 * atan(1.0));
  physics->H      = 5.0;
  physics->D      = 5.0;
  physics->Pm_avg = 0.9;
  physics->Pmax   = 2.1;
  physics->sigma  = 0.2;
  physics->lambda = 100.0 / physics->omegaB;

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
          if      (!strcmp(word, "omegaS")) {ferr=fscanf(in,"%lf",&physics->omegaS);if(ferr!=1)return(1);}
          else if (!strcmp(word, "omegaB")) {ferr=fscanf(in,"%lf",&physics->omegaB);if(ferr!=1)return(1);}
          else if (!strcmp(word, "H"     )) {ferr=fscanf(in,"%lf",&physics->H     );if(ferr!=1)return(1);}
          else if (!strcmp(word, "D"     )) {ferr=fscanf(in,"%lf",&physics->D     );if(ferr!=1)return(1);}
          else if (!strcmp(word, "Pm_avg")) {ferr=fscanf(in,"%lf",&physics->Pm_avg);if(ferr!=1)return(1);}
          else if (!strcmp(word, "Pmax"  )) {ferr=fscanf(in,"%lf",&physics->Pmax  );if(ferr!=1)return(1);}
          else if (!strcmp(word, "sigma" )) {ferr=fscanf(in,"%lf",&physics->sigma );if(ferr!=1)return(1);}
          else if (!strcmp(word, "lambda")) {ferr=fscanf(in,"%lf",&physics->lambda);if(ferr!=1)return(1);}
          else if (strcmp(word,"end")) {
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
  IERR MPIBroadcast_double(&physics->omegaS,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->omegaB,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->H     ,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->D     ,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->Pm_avg,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->Pmax  ,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->sigma ,1,0,&mpi->world);                        CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->lambda,1,0,&mpi->world);                        CHECKERR(ierr);
#endif

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in FPPowerSystem1BusInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = FPPowerSystem1BusComputeCFL;
  solver->ComputeDiffNumber  = FPPowerSystem1BusComputeDiffNumber;
  solver->FFunction          = FPPowerSystem1BusAdvection;
  solver->HFunction          = FPPowerSystem1BusDiffusion;
  solver->Upwind             = FPPowerSystem1BusUpwind;

  /* check that solver is using the correct diffusion formulation */
  if (strcmp(solver->spatial_type_par,_NC_2STAGE_)) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in FPPowerSystem1BusInitialize(): Parabolic term spatial discretization must be \"%s\"\n",_NC_2STAGE_);
    }
    return(1);
  }

  return(0);
}
