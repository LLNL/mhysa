#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

/* include header files for each physical model */
#include <linearadr.h>
#include <fokkerplanck.h>

int InitializePhysics(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0;

  if (!mpi->rank) printf("Initializing physics.\n");

  if (!strcmp(solver->model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {
    solver->physics = (LinearADR*) calloc (1,sizeof(LinearADR));
    ierr = LinearADRInitialize(solver,mpi); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_FOKKER_PLANCK_)) {
    solver->physics = (FokkerPlanck*) calloc (1,sizeof(FokkerPlanck));
    ierr = FokkerPlanckInitialize(solver,mpi); CHECKERR(ierr);
  } else {
    fprintf(stderr,"Error: %s is not a supported physical model.\n",solver->model);
    fprintf(stderr,"See header file \"physics.h\" for a list of supported models.\n");
    return(1);
  }

  return(0);
}
