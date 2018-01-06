/*! @file NavierStokes3DPreStep.c
    @brief Pre-step function for 3D Navier Stokes equations
    @author Debojyoti Ghosh
*/
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*! Pre-step function for the 3D Navier Stokes equations: This function
    is called at the beginning of each time step.
*/
int NavierStokes3DPreStep(
                            double  *u,   /*!< Solution vector */
                            void    *s,   /*!< Solver object of type #HyPar */
                            void    *m,   /*!< MPI object of type #MPIVariables */
                            double  waqt  /*!< Current simulation time */
                         )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;

  return(0);
}
