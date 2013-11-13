#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystemComputeCFL        (void*,void*,double,double);
double FPPowerSystemComputeDiffNumber (void*,void*,double,double);
int    FPPowerSystemAdvection         (double*,double*,int,void*,double);
int    FPPowerSystemDiffusion         (double*,double*,int,void*,double);
int    FPPowerSystemUpwind            (double*,double*,double*,double*,
                                       double*,double*,int,void*,double);
int    FPPowerSystemPostStep          (double*,void*,void*,double);
int    FPPowerSystemPrintStep         (void*,void*,double);

int FPPowerSystemInitialize(void *s,void *m)
{
  HyPar           *solver  = (HyPar*)         s;
  MPIVariables    *mpi     = (MPIVariables*)  m; 
  FPPowerSystem   *physics = (FPPowerSystem*) solver->physics;
  int             ferr;
  _DECLARE_IERR_;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in FPPowerSystemInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in FPPowerSystemInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
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
        if      (!strcmp(word, "omega_s")) {ferr=fscanf(in,"%lf",&physics->O_s);if(ferr!=1)return(1);}
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

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = FPPowerSystemComputeCFL;
  solver->ComputeDiffNumber  = FPPowerSystemComputeDiffNumber;
  solver->FFunction          = FPPowerSystemAdvection;
  solver->GFunction          = FPPowerSystemDiffusion;
  solver->Upwind             = FPPowerSystemUpwind;
  solver->PostStep           = FPPowerSystemPostStep;
  solver->PrintStep          = FPPowerSystemPrintStep;

  /* Calculate and print the PDF integral of the initial solution */
  IERR FPPowerSystemPostStep(solver->u,solver,mpi,0.0);        CHECKERR(ierr);
  if (!mpi->rank) IERR FPPowerSystemPrintStep(solver,mpi,0.0); CHECKERR(ierr);
  
  return(0);
}
