#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

double LinearADRComputeCFL        (void*,void*,double);
double LinearADRComputeDiffNumber (void*,void*,double);
int    LinearADRAdvection         (double*,double*,int,void*,double);
int    LinearADRDiffusion         (double*,double*,int,void*,double);
int    LinearADRReaction          ();
int    LinearADRUpwind            (double*,double*,double*,double*,int,void*);

int LinearADRInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  LinearADR     *physics = (LinearADR*)     solver->physics;
  int           ierr     = 0,i;

  physics->a = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));
  physics->d = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));

  /* default values are zero */
  ierr = ArraySetValue_double(physics->a,solver->ndims*solver->nvars,0.0); CHECKERR(ierr);
  ierr = ArraySetValue_double(physics->d,solver->ndims*solver->nvars,0.0); CHECKERR(ierr);

  /* reading physical model specific inputs - all processes */
  FILE *in;
  if (!mpi->rank) printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) {
    fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
    char word[_MAX_STRING_SIZE_];
    ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
        if (!strcmp(word, "advection")) {
          /* read advection coefficients */
          for (i=0; i<solver->ndims*solver->nvars; i++) ierr = fscanf(in,"%lf",&physics->a[i]);
          if (ierr != 1) return(1);
        } else if (!strcmp(word, "diffusion")) {
          /* read diffusion coefficients */
          for (i=0; i<solver->ndims*solver->nvars; i++) ierr = fscanf(in,"%lf",&physics->d[i]);
          if (ierr != 1) return(1);
        } else if (strcmp(word,"end")) {
          char useless[_MAX_STRING_SIZE_];
          ierr = fscanf(in,"%s",useless); if (ierr != 1) return(ierr);
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

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = LinearADRComputeCFL;
  solver->ComputeDiffNumber  = LinearADRComputeDiffNumber;
  solver->FFunction          = LinearADRAdvection;
  solver->GFunction          = LinearADRDiffusion;
  solver->SFunction          = LinearADRReaction;
  solver->Upwind             = LinearADRUpwind;

  return(0);
}
