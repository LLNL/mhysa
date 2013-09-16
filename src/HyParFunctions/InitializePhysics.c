#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>
#include <physics.h>

/* include header files for each physical model */
#include <advectiondiffusionreaction.h>

int InitializePhysics(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0;

  /* Initialize all physics-model related function pointers to NULL */
  solver->ComputeCFL          = NULL;
  solver->ComputeDiffNumber   = NULL;
  solver->HyperbolicFunction  = NULL;
  solver->ParabolicFunction   = NULL;
  solver->SourceFunction      = NULL; 

  if (!mpi->rank) printf("Initializing physics.\n");

  if (!strcmp(solver->model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {
    ierr = LinearADRInitialize(solver,mpi); CHECKERR(ierr);
    solver->ComputeCFL          = LinearADRComputeCFL;
    solver->ComputeDiffNumber   = LinearADRComputeDiffNumber;
    solver->HyperbolicFunction  = LinearADRAdvection;
    solver->ParabolicFunction   = LinearADRDiffusion;
    solver->SourceFunction      = LinearADRReaction;
  } else {
    fprintf(stderr,"Error: %s is not a supported physical model.\n",solver->model);
    fprintf(stderr,"See header file \"physics.h\" for a list of supported models.\n");
    return(1);
  }

  return(0);
}
