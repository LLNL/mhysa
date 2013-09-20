#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

//double FPDoubleWellComputeCFL        (void*,void*,double);
//double FPDoubleWellComputeDiffNumber (void*,void*,double);
//int    FPDoubleWellAdvection         (double*,double*,int,void*);
//int    FPDoubleWellDiffusion         (double*,double*,int,void*);
//int    FPDoubleWellUpwind            (double*,double*,double*,double*,int,void*);

int FPDoubleWellInitialize(void *s,void *m)
{
  HyPar         *solver  = (HyPar*)         s;
//  MPIVariables  *mpi     = (MPIVariables*)  m; 
//  FPDoubleWell  *physics = (FPDoubleWell*)  solver->physics;
//  int           ierr     = 0,i;

  if (solver->nvars != 1) {
    fprintf(stderr,"Error in FPDoubleWellInitializeO(): nvars has to be 1.\n");
    return(1);
  }

//  physics->drift = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));
//  physics->diff  = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));

  /* default values are zero */
//  ierr = ArraySetValue_double(physics->a,solver->ndims*solver->nvars,0.0); CHECKERR(ierr);
//  ierr = ArraySetValue_double(physics->d,solver->ndims*solver->nvars,0.0); CHECKERR(ierr);

  /* initializing physical model-specific functions */
//  solver->ComputeCFL         = FPDoubleWellComputeCFL;
//  solver->ComputeDiffNumber  = FPDoubleWellComputeDiffNumber;
//  solver->FFunction          = FPDoubleWellAdvection;
//  solver->GFunction          = FPDoubleWellDiffusion;
//  solver->SFunction          = NULL;
//  solver->Upwind             = FPDoubleWellUpwind;

  return(0);
}
