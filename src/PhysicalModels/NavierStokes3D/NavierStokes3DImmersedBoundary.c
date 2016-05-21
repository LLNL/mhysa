/*! @file NavierStokes3DImmersedBoundary.c
    @brief Immersed boundary treatment for 3D Navier-Stokes equations
    @author Debojyoti Ghosh
*/

#include <arrayfunctions.h>
#include <immersedboundaries.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*!
  Apply no-slip wall boundary conditions on the immersed boundary
  points (grid points within the immersed body that are within 
  stencil-width distance of interior points, i.e., points in the 
  interior of the computational domain).
*/
int NavierStokes3DImmersedBoundary(
                                    void    *s, /*!< Solver object of type #HyPar */
                                    double  *u, /*!< Array with the solution vector */ 
                                    double  t   /*!< Current simulation time */
                                  )
{
  HyPar             *solver   = (HyPar*)   s;
  ImmersedBoundary  *IB       = (ImmersedBoundary*) solver->ib;
  IBNode            *boundary = IB->boundary;
  NavierStokes3D    *param    = (NavierStokes3D*) solver->physics;
  static double     v[_MODEL_NVARS_];
  int               n, j, k, nb = IB->n_boundary_nodes;

  if (!solver->flag_ib) return(0);

  double inv_gamma_m1 = 1.0 / (param->gamma - 1.0);

  for (n=0; n<nb; n++) {

    int     node_index = boundary[n].p;
    double  *alpha = &(boundary[n].interp_coeffs[0]);
    int     *nodes = &(boundary[n].interp_nodes[0]);
    double  factor = boundary[n].surface_distance / boundary[n].interp_node_distance;

    _ArraySetValue_(v,_MODEL_NVARS_,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<_MODEL_NVARS_; k++) {
        v[k] += ( alpha[j] * u[_MODEL_NVARS_*nodes[j]+k] );
      }
    }

    double rho, uvel, vvel, wvel, energy, pressure;
    _NavierStokes3DGetFlowVar_(v,rho,uvel,vvel,wvel,energy,pressure,param);

    double rho_ib, uvel_ib, vvel_ib, wvel_ib, energy_ib, pressure_ib;
    rho_ib = rho;
    pressure_ib = pressure;
    uvel_ib = -uvel * factor;
    vvel_ib = -vvel * factor;
    wvel_ib = -wvel * factor;
    energy_ib = inv_gamma_m1*pressure_ib + 0.5*rho_ib*(uvel_ib*uvel_ib+vvel_ib*vvel_ib+wvel_ib*wvel_ib);

    u[_MODEL_NVARS_*node_index+0] = rho_ib;
    u[_MODEL_NVARS_*node_index+1] = rho_ib * uvel_ib;
    u[_MODEL_NVARS_*node_index+2] = rho_ib * vvel_ib;
    u[_MODEL_NVARS_*node_index+3] = rho_ib * wvel_ib;
    u[_MODEL_NVARS_*node_index+4] = energy_ib;
  }
  return(0);
}
