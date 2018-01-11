/*! @file InitializePhysics.c
    @author Debojyoti Ghosh
    @brief Initialize the physical model
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <bandedmatrix.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* include header files for each physical model */
#include <physicalmodels/euler1d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Initialize the physical model for a simulation: Depending on the 
    physical model specified, this function calls the initialization
    function for that physical model. The latter is responsible for
    setting all the physics-specific functions that are required
    by the model.
*/
int InitializePhysics(
                        void *s, /*!< Solver object of type #HyPar */
                        void *m  /*!< MPI object of type #MPIVariables */
                     )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  if (!mpi->rank) printf("Initializing physics. Model = \"%s\"\n",solver->model);

  /* Initialize physics-specific functions to NULL */
  solver->ComputeCFL            = NULL;
  solver->ComputeDiffNumber     = NULL;
  solver->FFunction             = NULL;
  solver->dFFunction            = NULL;
  solver->FdFFunction           = NULL;
  solver->GFunction             = NULL;
  solver->HFunction             = NULL;
  solver->SFunction             = NULL;
  solver->UFunction             = NULL;
  solver->JFunction             = NULL;
  solver->Upwind                = NULL;
  solver->UpwinddF              = NULL;
  solver->UpwindFdF             = NULL;
  solver->PreStage              = NULL;
  solver->PostStage             = NULL;
  solver->PreStep               = NULL;
  solver->PostStep              = NULL;
  solver->PrintStep             = NULL;
  solver->PhysicsOutput         = NULL;
  solver->AveragingFunction     = NULL;
  solver->GetLeftEigenvectors   = NULL;
  solver->GetRightEigenvectors  = NULL;
  solver->IBFunction            = NULL;

  if (!strcmp(solver->model,_EULER_1D_)) {

    solver->physics = (Euler1D*) calloc (1,sizeof(Euler1D));
    IERR Euler1DInitialize(solver,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->model,_NAVIER_STOKES_3D_)) {

    solver->physics = (NavierStokes3D*) calloc (1,sizeof(NavierStokes3D));
    IERR NavierStokes3DInitialize(solver,mpi); CHECKERR(ierr);

  } else {

    fprintf(stderr,"Error: %s is not a supported physical model.\n",solver->model);
    return(1);

  }

  if (solver->Upwind == NULL) {
    return(1);
  }

  /* some checks */
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

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if ((!solver->dFFunction) || (!solver->UpwinddF)) {
      if (!mpi->rank) {
        fprintf(stderr,"Error: Splitting of hyperbolic flux requires a dFFunction ");
        fprintf(stderr,"and its upwinding function UpwinddF.\n");
        fprintf(stderr,"Error: f(u) = [f(u) - df(u)] + df(u).\n");
        fprintf(stderr,"Error: dFFunction or UpwinddF (or both) is (are) NULL.\n");
      }
      return(1);
    }
    if (solver->FdFFunction && solver->UpwindFdF) solver->flag_fdf_specified = 1;
    else                                          solver->flag_fdf_specified = 0;
  }

  if ((solver->IBFunction == NULL) && (solver->flag_ib)) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in InitializePhysics(): Physical model %s does not yet have an immersed boundary treatment.\n",
              solver->model);
    }
    return(1);
  }

  return(0);
}
