#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* include header files for each physical model */
#include <physicalmodels/linearadr.h>
#include <physicalmodels/fpdoublewell.h>
#include <physicalmodels/fppowersystem.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <physicalmodels/euler1d.h>
#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes2d.h>
#include <physicalmodels/navierstokes3d.h>
#include <physicalmodels/numa3d.h>

int InitializePhysics(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  if (!mpi->rank) printf("Initializing physics. Model = \"%s\"\n",solver->model);

  /* Initialize physics-specific functions to NULL */
  solver->ComputeCFL            = NULL;
  solver->ComputeDiffNumber     = NULL;
  solver->FFunction             = NULL;
  solver->GFunction             = NULL;
  solver->HFunction             = NULL;
  solver->SFunction             = NULL;
  solver->Upwind                = NULL;
  solver->PreStage              = NULL;
  solver->PostStage             = NULL;
  solver->PreStep               = NULL;
  solver->PostStep              = NULL;
  solver->PrintStep             = NULL;
  solver->AveragingFunction     = NULL;
  solver->GetLeftEigenvectors   = NULL;
  solver->GetRightEigenvectors  = NULL;

  if (!strcmp(solver->model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {

    solver->physics = (LinearADR*) calloc (1,sizeof(LinearADR));
    IERR LinearADRInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_FP_DOUBLE_WELL_)) {

    solver->physics = (FPDoubleWell*) calloc (1,sizeof(FPDoubleWell));
    IERR FPDoubleWellInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_FP_POWER_SYSTEM_)) {

    solver->physics = (FPPowerSystem*) calloc (1,sizeof(FPPowerSystem));
    IERR FPPowerSystemInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_FP_POWER_SYSTEM_3BUS_)) {

    solver->physics = (FPPowerSystem3Bus*) calloc (1,sizeof(FPPowerSystem3Bus));
    IERR FPPowerSystem3BusInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_EULER_1D_)) {

    solver->physics = (Euler1D*) calloc (1,sizeof(Euler1D));
    IERR Euler1DInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_EULER_2D_)) {

    solver->physics = (Euler2D*) calloc (1,sizeof(Euler2D));
    IERR Euler2DInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_NAVIER_STOKES_2D_)) {

    solver->physics = (NavierStokes2D*) calloc (1,sizeof(NavierStokes2D));
    IERR NavierStokes2DInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_NAVIER_STOKES_3D_)) {

    solver->physics = (NavierStokes3D*) calloc (1,sizeof(NavierStokes3D));
    IERR NavierStokes3DInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_NUMA3D_)) {

    solver->physics = (Numa3D*) calloc (1,sizeof(Numa3D));
    IERR Numa3DInitialize(solver,mpi); CHECKERR(ierr);

  } else {

    fprintf(stderr,"Error: %s is not a supported physical model.\n",solver->model);
    return(1);

  }

  if ( ( (solver->GetLeftEigenvectors == NULL) || (solver->GetRightEigenvectors == NULL) )
      && (!strcmp(solver->interp_type,_CHARACTERISTIC_)) && (solver->nvars > 1) ) {
    if (!mpi->rank) {
      fprintf(stderr,"Error: Interpolation type is defined as characteristic ");
      fprintf(stderr,"but physics initializations returned NULL pointers for ");
      fprintf(stderr,"Get(Left,Right)Eigenvectors needed for characteristic-based ");
      fprintf(stderr,"reconstruction.\n");
    }
    return(1);
  }


  return(0);
}
