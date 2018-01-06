/*! @file BCSubsonicOutflow.c
    @author Debojyoti Ghosh
    @brief Subsonic outflow boundary conditions (specific to Euler/Navier-Stokes systems)
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/navierstokes3d.h>

/*! Applies the subsonic outflow boundary condition: The pressure at the physical boundary
    ghost points is specified, while the density and velocity are extrapolated from the 
    interior. This boundary condition is specific to the 3D Navier-Stokes system 
    (#NavierStokes3D).
*/
int BCSubsonicOutflowU(
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

  int dim   = boundary->dim;
  int face  = boundary->face;
  int k;

  if (ndims == 3) {

    NavierStokes3D *physics = (NavierStokes3D*) (*(boundary->physics));
    double gamma = physics->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);
    int ns = physics->n_species;
    int nv = physics->n_vibeng;

    if (boundary->on_this_proc) {

      int bounds[ndims], indexb[ndims], indexi[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);

      int done = 0;
      while (!done) {

        int p1, p2;
        _ArrayCopy1D_(indexb,indexi,ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
        if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
        
        /* flow variables in the interior */
        double rho_s[ns], rho_t, uvel, vvel, wvel, E, E_v[nv], pressure, T;
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho_s,rho_t,uvel,vvel,wvel,E,E_v,pressure,T,physics);

        double rho_s_gpt[ns], rho_t_gpt, uvel_gpt, vvel_gpt, wvel_gpt, E_gpt, E_v_gpt[nv],pressure_gpt;
        for (k = 0; k < ns; k++) rho_s_gpt[k] = rho_s[k];
        rho_t_gpt = rho_t;
        pressure_gpt = pressure; /* useless statement to avoid compiler warning */
        pressure_gpt = boundary->FlowPressure;
        uvel_gpt = uvel;
        vvel_gpt = vvel;
        wvel_gpt = wvel;
        E_gpt = inv_gamma_m1*pressure_gpt/rho_t_gpt
                + 0.5 * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);
        for (k = 0; k < nv; k++) E_v_gpt[k] = E_v[k];

        _NavierStokes3DSetFlowVar_( (phi+nvars*p1), rho_s_gpt, rho_t_gpt,
                                    uvel_gpt, vvel_gpt, wvel_gpt,
                                    E_gpt, E_v_gpt, pressure_gpt, physics );

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  
  return(0);
}
