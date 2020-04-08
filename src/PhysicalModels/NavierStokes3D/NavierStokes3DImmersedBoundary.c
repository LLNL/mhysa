/*! @file NavierStokes3DImmersedBoundary.c
    @brief Immersed boundary treatment for 3D Navier-Stokes equations
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <immersedboundaries.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Apply no-slip adiabatic wall boundary conditions on the immersed boundary
    points (grid points within the immersed body that are within 
    stencil-width distance of interior points, i.e., points in the 
    interior of the computational domain). */
int NavierStokes3DIBAdiabatic(void    *s, /*!< Solver object of type #HyPar */
                              void    *m, /*!< Solver object of type #HyPar */
                              double  *u, /*!< Array with the solution vector */ 
                              double  t   /*!< Current simulation time */
                             )
{
  HyPar             *solver   = (HyPar*)   s;
  MPIVariables      *mpi      = (MPIVariables*) m;
  ImmersedBoundary  *IB       = (ImmersedBoundary*) solver->ib;
  IBNode            *boundary = IB->boundary;
  NavierStokes3D    *param    = (NavierStokes3D*) solver->physics;
  
  int nvars = solver->nvars;
  int ns    = param->n_species;
  int nv    = param->n_vibeng;

  double  v[nvars];
  int     n, j, k, nb = IB->n_boundary_nodes;

  if (!solver->flag_ib) return(0);

  /* Ideally, this shouldn't be here - But this function is called everywhere
     (through ApplyIBConditions()) *before* MPIExchangeBoundariesnD is called! */
  MPIExchangeBoundariesnD(  _MODEL_NDIMS_,
                            nvars,
                            solver->dim_local,
                            solver->ghosts,
                            mpi,
                            u );

  double inv_gamma_m1 = 1.0 / (param->gamma - 1.0);

  for (n=0; n<nb; n++) {

    int     node_index = boundary[n].p;
    double  *alpha = &(boundary[n].interp_coeffs[0]);
    int     *nodes = &(boundary[n].interp_nodes[0]);
    double  factor = boundary[n].surface_distance / boundary[n].interp_node_distance;

    _ArraySetValue_(v,nvars,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<nvars; k++) {
        v[k] += ( alpha[j] * u[nvars*nodes[j]+k] );
      }
    }

    double rho_s[ns], rho_t, uvel, vvel, wvel, E, E_v[nv], pressure, T;
    _NavierStokes3DGetFlowVar_(v,rho_s,rho_t,uvel,vvel,wvel,E,E_v,pressure,T,param);

    double rho_s_ib[ns], rho_t_ib, uvel_ib, vvel_ib, wvel_ib, E_ib, E_v_ib[nv], pressure_ib;
    for (k = 0; k < ns; k++)  rho_s_ib[k] = rho_s[k];
    _NavierStokes3DTotalDensity_(rho_t_ib,rho_s_ib,ns)
    pressure_ib = pressure;
    uvel_ib = -uvel * factor;
    vvel_ib = -vvel * factor;
    wvel_ib = -wvel * factor;
    E_ib = inv_gamma_m1*pressure_ib/rho_t_ib 
            + 0.5*(uvel_ib*uvel_ib+vvel_ib*vvel_ib+wvel_ib*wvel_ib);

    for (k = 0; k < nv; k++)  E_v_ib[k] = E_v[k];
    _NavierStokes3DSetFlowVar_( (u+nvars*node_index),rho_s_ib,rho_t_ib,
                                uvel_ib,vvel_ib,wvel_ib,
                                E_ib,E_v_ib,pressure_ib,param );
  }
  return(0);
}

/*! Apply no-slip isothermal wall boundary conditions on the immersed boundary
    points (grid points within the immersed body that are within 
    stencil-width distance of interior points, i.e., points in the 
    interior of the computational domain). */
int NavierStokes3DIBIsothermal( void    *s, /*!< Solver object of type #HyPar */
                                void    *m, /*!< Solver object of type #HyPar */
                                double  *u, /*!< Array with the solution vector */ 
                                double  t   /*!< Current simulation time */
                              )
{
  HyPar             *solver   = (HyPar*)   s;
  MPIVariables      *mpi      = (MPIVariables*) m;
  ImmersedBoundary  *IB       = (ImmersedBoundary*) solver->ib;
  IBNode            *boundary = IB->boundary;
  NavierStokes3D    *param    = (NavierStokes3D*) solver->physics;
  
  int nvars = solver->nvars;
  int ns    = param->n_species;
  int nv    = param->n_vibeng;

  double v[nvars];
  int    n, j, k, nb = IB->n_boundary_nodes;

  if (!solver->flag_ib) return(0);

  /* Ideally, this shouldn't be here - But this function is called everywhere
     (through ApplyIBConditions()) *before* MPIExchangeBoundariesnD is called! */
  MPIExchangeBoundariesnD(  _MODEL_NDIMS_,
                            nvars,
                            solver->dim_local,
                            solver->ghosts,
                            mpi,
                            u );

  double inv_gamma_m1 = 1.0 / (param->gamma - 1.0);

  for (n=0; n<nb; n++) {

    int     node_index = boundary[n].p;
    double  *alpha = &(boundary[n].interp_coeffs[0]);
    int     *nodes = &(boundary[n].interp_nodes[0]);
    double  factor = boundary[n].surface_distance / boundary[n].interp_node_distance;

    _ArraySetValue_(v,nvars,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<nvars; k++) {
        v[k] += ( alpha[j] * u[nvars*nodes[j]+k] );
      }
    }

    double rho_s[ns], rho_t, uvel, vvel, wvel, E, E_v[nv], pressure, T;
    _NavierStokes3DGetFlowVar_(v,rho_s,rho_t,uvel,vvel,wvel,E,E_v,pressure,T,param);
    double mfrac[ns];
    for (k = 0; k < ns; k++) mfrac[k] = rho_s[k]/rho_t;

    double rho_s_ib[ns], rho_t_ib, uvel_ib, vvel_ib, wvel_ib, E_ib, E_v_ib[nv], pressure_ib, T_ib;
    pressure_ib = pressure;
    T_ib = (1.0+factor)*param->T_ib_wall - factor * T;
    if (T_ib < _MACHINE_ZERO_) T_ib = param->T_ib_wall;
    rho_t_ib = pressure_ib / T_ib;
    for (k = 0; k < ns; k++) rho_s_ib[k] = mfrac[k] * rho_t_ib;
    uvel_ib = -uvel * factor;
    vvel_ib = -vvel * factor;
    wvel_ib = -wvel * factor;
    E_ib = inv_gamma_m1*pressure_ib/rho_t_ib 
            + 0.5*(uvel_ib*uvel_ib+vvel_ib*vvel_ib+wvel_ib*wvel_ib);

    for (k = 0; k < nv; k++)  E_v_ib[k] = E_v[k];
    _NavierStokes3DSetFlowVar_( (u+nvars*node_index),rho_s_ib,rho_t_ib,
                                uvel_ib,vvel_ib,wvel_ib,
                                E_ib,E_v_ib,pressure_ib,param );
  }

  return(0);
}
