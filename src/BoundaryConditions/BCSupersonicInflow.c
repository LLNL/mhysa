/*! @file BCSupersonicInflow.c
    @author Debojyoti Ghosh
    @brief Supersonic inflow boundary conditions (specific to Euler/Navier-Stokes systems)
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/navierstokes3d.h>

/*! Applies the supersonic (steady) inflow boundary condition: All the flow variables
    (density, pressure, velocity) are specified at the physical boundary ghost points,
    since it is supersonic inflow. This boundary condition is specific to the 3D 
    Navier-Stokes systems (#NavierStokes3D).
    \n\n
    Note: the Dirichlet boundary condition (#_DIRICHLET_) could be used as well for
    supersonic inflow; however the specified Dirichlet state should be in terms of the 
    conserved variables, while the specified supersonic inflow state here is in terms of
    the flow variables.
*/
int BCSupersonicInflowU(
                        void    *b,     /*!< Boundary object of type #DomainBoundary */
                        void    *m,     /*!< MPI object of type #MPIVariables */
                        int     ndims,  /*!< Number of spatial dimensions */
                        int     nvars,  /*!< Number of variables/DoFs per grid point */
                        int     *size,  /*!< Integer array with the number of grid points in each spatial dimension */
                        int     ghosts, /*!< Number of ghost points */
                        double  *phi,   /*!< The solution array on which to apply the boundary condition */
                        double  waqt    /*!< Current solution time */
                       )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  if (ndims == 3) {

    NavierStokes3D* physics = (NavierStokes3D*) (*(boundary->physics));
    double gamma = physics->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->on_this_proc) {
      int bounds[ndims], indexb[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);
      int done = 0;
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt      = boundary->FlowDensity;
        pressure_gpt = boundary->FlowPressure;
        uvel_gpt     = boundary->FlowVelocity[0];
        vvel_gpt     = boundary->FlowVelocity[1];
        wvel_gpt     = boundary->FlowVelocity[2];
        energy_gpt   = inv_gamma_m1*pressure_gpt
                       + 0.5 * rho_gpt 
                       * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);

        phi[nvars*p1+0] = rho_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt;
        phi[nvars*p1+3] = rho_gpt * wvel_gpt;
        phi[nvars*p1+4] = energy_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }

  return(0);
}
