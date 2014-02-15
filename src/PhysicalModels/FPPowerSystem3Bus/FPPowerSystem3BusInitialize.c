#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystem3BusComputeCFL        (void*,void*,double,double);
double FPPowerSystem3BusComputeDiffNumber (void*,void*,double,double);
int    FPPowerSystem3BusAdvection         (double*,double*,int,void*,double);
int    FPPowerSystem3BusDiffusion         (double*,double*,int,void*,double);
int    FPPowerSystem3BusUpwind            (double*,double*,double*,double*,
                                           double*,double*,int,void*,double);
int    FPPowerSystem3BusPostStep          (double*,void*,void*,double);
int    FPPowerSystem3BusPrintStep         (void*,void*,double);

int FPPowerSystem3BusInitialize(void *s,void *m)
{
  HyPar               *solver  = (HyPar*)             s;
  MPIVariables        *mpi     = (MPIVariables*)      m; 
//  FPPowerSystem3Bus   *physics = (FPPowerSystem3Bus*) solver->physics;
  int                 ferr;
  _DECLARE_IERR_;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values of model parameters */

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
      }
	  } else {
    	fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
  }
  fclose(in);

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = FPPowerSystem3BusComputeCFL;
  solver->ComputeDiffNumber  = FPPowerSystem3BusComputeDiffNumber;
  solver->FFunction          = FPPowerSystem3BusAdvection;
  solver->GFunction          = FPPowerSystem3BusDiffusion;
  solver->Upwind             = FPPowerSystem3BusUpwind;
  solver->PostStep           = FPPowerSystem3BusPostStep;
  solver->PrintStep          = FPPowerSystem3BusPrintStep;

  /* Calculate and print the PDF integral of the initial solution */
  IERR FPPowerSystem3BusPostStep(solver->u,solver,mpi,0.0);        CHECKERR(ierr);
  if (!mpi->rank) IERR FPPowerSystem3BusPrintStep(solver,mpi,0.0); CHECKERR(ierr);
  
  return(0);
}
