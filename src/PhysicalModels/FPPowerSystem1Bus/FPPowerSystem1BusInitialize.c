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
int    FPPowerSystem1BusDiffusion         (double*,double*,int,void*,double);
int    FPPowerSystem1BusUpwind            (double*,double*,double*,double*,
                                       double*,double*,int,void*,double);
int    FPPowerSystem1BusPostStep          (double*,void*,void*,double);
int    FPPowerSystem1BusPrintStep         (void*,void*,double);

int FPPowerSystem1BusInitialize(void *s,void *m)
{
  HyPar           *solver  = (HyPar*)         s;
  MPIVariables    *mpi     = (MPIVariables*)  m; 
  FPPowerSystem1Bus   *physics = (FPPowerSystem1Bus*) solver->physics;
  int             ferr;
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
  physics->H   = 5.0;
  physics->D   = 5.0;
  physics->E   = 1.1358;
  physics->V   = 1.0;
  physics->g1  = 0.545;
  physics->g2  = 0.745;
  physics->Pm  = 0.9;
  physics->l   = 0.1;
  physics->q   = 1.0;
  physics->O_s = 120*(4.0*atan(1.0));
  physics->tf  = 0.1;
  physics->tcl = 0.2;

  /* reading physical model specific inputs - all processes */
  FILE *in;
  if (!mpi->rank) printf("Reading physical model inputs from file \"physics.inp\".\n");
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
        if      (!strcmp(word, "inertia")) {ferr=fscanf(in,"%lf",&physics->H  );if(ferr!=1)return(1);}
        else if (!strcmp(word, "omega_s")) {ferr=fscanf(in,"%lf",&physics->O_s);if(ferr!=1)return(1);}
        else if (!strcmp(word, "E"      )) {ferr=fscanf(in,"%lf",&physics->E  );if(ferr!=1)return(1);}
        else if (!strcmp(word, "V"      )) {ferr=fscanf(in,"%lf",&physics->V  );if(ferr!=1)return(1);}
        else if (!strcmp(word, "g1"     )) {ferr=fscanf(in,"%lf",&physics->g1 );if(ferr!=1)return(1);}
        else if (!strcmp(word, "g2"     )) {ferr=fscanf(in,"%lf",&physics->g2 );if(ferr!=1)return(1);}
        else if (!strcmp(word, "D"      )) {ferr=fscanf(in,"%lf",&physics->D  );if(ferr!=1)return(1);}
        else if (!strcmp(word, "PM_min" )) {ferr=fscanf(in,"%lf",&physics->Pm );if(ferr!=1)return(1);}
        else if (!strcmp(word, "lambda" )) {ferr=fscanf(in,"%lf",&physics->l  );if(ferr!=1)return(1);}
        else if (!strcmp(word, "q"      )) {ferr=fscanf(in,"%lf",&physics->q  );if(ferr!=1)return(1);}
        else if (!strcmp(word, "tf"     )) {ferr=fscanf(in,"%lf",&physics->tf );if(ferr!=1)return(1);}
        else if (!strcmp(word, "tcl"    )) {ferr=fscanf(in,"%lf",&physics->tcl);if(ferr!=1)return(1);}
      }
	  } else {
    	fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
  }
  fclose(in);

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
  solver->GFunction          = FPPowerSystem1BusDiffusion;
  solver->Upwind             = FPPowerSystem1BusUpwind;
  solver->PostStep           = FPPowerSystem1BusPostStep;
  solver->PrintStep          = FPPowerSystem1BusPrintStep;

  /* Calculate and print the PDF integral of the initial solution */
  IERR FPPowerSystem1BusPostStep(solver->u,solver,mpi,0.0);        CHECKERR(ierr);
  if (!mpi->rank) IERR FPPowerSystem1BusPrintStep(solver,mpi,0.0); CHECKERR(ierr);
  
  return(0);
}
