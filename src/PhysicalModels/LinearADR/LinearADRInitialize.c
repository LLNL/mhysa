#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

double LinearADRComputeCFL        (void*,void*,double,double);
double LinearADRComputeDiffNumber (void*,void*,double,double);
int    LinearADRAdvection         (double*,double*,int,void*,double);
int    LinearADRDiffusion         (double*,double*,int,void*,double);
int    LinearADRReaction          ();
int    LinearADRUpwind            (double*,double*,double*,double*,
                                   double*,double*,int,void*,double);

int LinearADRInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  LinearADR     *physics = (LinearADR*)     solver->physics;
  int           i,ferr;

  physics->a = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));
  physics->d = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));

  /* default values are zero */
  _ArraySetValue_(physics->a,solver->ndims*solver->nvars,0.0);
  _ArraySetValue_(physics->d,solver->ndims*solver->nvars,0.0);

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
          if (!strcmp(word, "advection")) {
            /* read advection coefficients */
            for (i=0; i<solver->ndims*solver->nvars; i++) ferr = fscanf(in,"%lf",&physics->a[i]);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "diffusion")) {
            /* read diffusion coefficients */
            for (i=0; i<solver->ndims*solver->nvars; i++) ferr = fscanf(in,"%lf",&physics->d[i]);
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
  IERR MPIBroadcast_double(physics->a,solver->ndims*solver->nvars,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double(physics->d,solver->ndims*solver->nvars,0,&mpi->world); CHECKERR(ierr);
#endif

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in LinearADRInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = LinearADRComputeCFL;
  solver->ComputeDiffNumber  = LinearADRComputeDiffNumber;
  solver->FFunction          = LinearADRAdvection;
  solver->GFunction          = LinearADRDiffusion;
  solver->SFunction          = LinearADRReaction;
  solver->Upwind             = LinearADRUpwind;

  return(0);
}
