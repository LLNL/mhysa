/*! @file NavierStokes2DPostStage.c
    @brief Pre-step function for 2D Navier Stokes equations
    @author Debojyoti Ghosh
*/
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*! Post-stage (multi-stage time integration)  function for the 2D Navier Stokes equations: 
    This function is called after each stage in a multi-stage time integration.
    + The stage solution is copied into #NavierStokes2D::solution for a linearized upwinding
      of partitioned flux terms.
*/
int NavierStokes2DPostStage(
                            double  *u,   /*!< Solution vector */
                            void    *s,   /*!< Solver object of type #HyPar */
                            void    *m,   /*!< MPI object of type #MPIVariables */
                            double  waqt  /*!< Current simulation time */
                         )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes2D    *param  = (NavierStokes2D*) solver->physics;

  /* copy the solution to act as a reference for linearization */
  _ArrayCopy1D_(u,param->solution,(solver->npoints_local_wghosts*_MODEL_NVARS_));

  return(0);
}
