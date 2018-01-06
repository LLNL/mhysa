/*! @file NavierStokes3DSource.c
    @author Debojyoti Ghosh
    @brief Compute the source term for the 3D Navier Stokes system
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the species density change rates, given species densities and temperature */
static void ComputeRates(
                          double  *rho_dot, /*!< Computed density change rates 
                                                 (array of size #NavierStokes3D::n_species */
                          double  *rho,     /*!< Species density (array of size #NavierStokes3D::n_species */
                          double  T,        /*!< Temperature */
                          void    *ctxt     /*!< Context object of type #NavierStokes3D */
                        )
{
  NavierStokes3D *param = (NavierStokes3D*) ctxt;
  int ns = param->n_species;

  /*
   * This is where the rate calculation will go in
   * rho_dot[0], rho_dot[1], ..., rho_dot[ns-1]
   * should be calculated from
   * rho[0], rho[1], ..., rho[ns-1], T
  */

  return;
}

/*! Computes the source term for the 3D Navier-Stokes equations
*/
int NavierStokes3DSource(
                          double  *source,  /*!< Array to hold the computed source */
                          double  *u,       /*!< Solution vector array */
                          void    *s,       /*!< Solver object of type #HyPar */
                          void    *m,       /*!< MPI object of type #MPIVariables */
                          double  t         /*!< Current simulation time */
                        )
{
  HyPar           *solver = (HyPar* ) s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ghosts  = solver->ghosts;
  int *dim    = solver->dim_local;
  int ns      = param->n_species;
  int nv      = param->n_vibeng;

  int     index[ndims], done;
  double  rho_s[ns], rho_t, vx, vy, vz, E, E_v[nv], P, T;
  double  rho_dot_s[ns], E_v_dot[nv];

  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    _NavierStokes3DGetFlowVar_((u+nvars*p),rho_s,rho_t,vx,vy,vz,E,E_v,P,T,param);
    _ArraySetValue_(rho_dot_s,ns,0.0);
    _ArraySetValue_(E_v_dot,nv,0.0);

    ComputeRates(rho_dot_s, rho_s, T, param);

    _NavierStokes3DSetSource_((source+nvars*p),rho_dot_s,E_v_dot,param);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(0);
}
