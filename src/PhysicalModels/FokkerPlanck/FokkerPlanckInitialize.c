#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <fokkerplanck.h>
#include <mpivars.h>
#include <hypar.h>

//double FokkerPlanckComputeCFL        (void*,void*,double);
//double FokkerPlanckComputeDiffNumber (void*,void*,double);
//int    FokkerPlanckAdvection         (double*,double*,int,void*);
//int    FokkerPlanckDiffusion         (double*,double*,int,void*);
//int    FokkerPlanckUpwind            (double*,double*,double*,double*,int,void*);

int FokkerPlanckInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
//  MPIVariables  *mpi     = (MPIVariables*)  m; 
//  FokkerPlanck  *physics = (FokkerPlanck*)  solver->physics;
//  int           ierr     = 0,i;

  if (solver->nvars != 1) {
    fprintf(stderr,"Error in FokkerPlanckInitializeO(): nvars has to be 1.\n");
    return(1);
  }

//  physics->drift = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));
//  physics->diff  = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));

  /* default values are zero */
//  ierr = ArraySetValue_double(physics->a,solver->ndims*solver->nvars,0.0); CHECKERR(ierr);
//  ierr = ArraySetValue_double(physics->d,solver->ndims*solver->nvars,0.0); CHECKERR(ierr);

  /* initializing physical model-specific functions */
//  solver->ComputeCFL         = FokkerPlanckComputeCFL;
//  solver->ComputeDiffNumber  = FokkerPlanckComputeDiffNumber;
//  solver->FFunction          = FokkerPlanckAdvection;
//  solver->GFunction          = FokkerPlanckDiffusion;
//  solver->SFunction          = NULL;
//  solver->Upwind             = FokkerPlanckUpwind;

  return(0);
}
