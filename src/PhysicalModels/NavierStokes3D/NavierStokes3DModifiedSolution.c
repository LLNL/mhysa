/*! @file NavierStokes3DModifiedSolution.c
    @author Debojyoti Ghosh
    @brief Compute the modified solution for the well-balanced treatment of gravitational source terms.
*/
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    This function computes the modified solution for the well-balanced treatment of the 
    gravitational source terms. The modified solution vector is given by
    \f{equation}{
      {\bf u}^* = \left[\begin{array}{c} \rho \varrho^{-1}\left(x,y\right) \\ \rho u \varrho^{-1}\left(x,y\right) \\ \rho v \varrho^{-1}\left(x,y\right) \\ \rho w \varrho^{-1}\left(x,y\right) \\ e^* \end{array}\right]
    \f}
    where 
    \f{equation}{
      e^* = \frac{p \varphi^{-1}\left(x,y\right)}{\gamma-1} + \frac{1}{2}\rho \varrho^{-1}\left(x,y\right) \left(u^2+v^2+w^2\right)
    \f}
    \f$\varrho\f$ and \f$\varphi\f$ are computed in #NavierStokes3DGravityField(). For flows without gravity, \f${\bf u}^* = {\bf u}\f$.
    
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.
*/
int NavierStokes3DModifiedSolution(
                                    double  *uC,  /*!< Array to hold the computed modified solution */
                                    double  *u,   /*!< Solution vector array */
                                    int     d,    /*!< spatial dimension (not used) */
                                    void    *s,   /*!< Solver object of type #HyPar */
                                    void    *m,   /*!< MPI object of time #MPIVariables */
                                    double waqt   /*!< Current simulation time */
                                  )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  double inv_gamma_m1 = 1.0 / (param->gamma-1.0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, uvel, vvel, wvel, E, P;
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,wvel,E,P,param);
    uC[_MODEL_NVARS_*p+0] = u[_MODEL_NVARS_*p+0] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+1] = u[_MODEL_NVARS_*p+1] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+2] = u[_MODEL_NVARS_*p+2] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+3] = u[_MODEL_NVARS_*p+3] * param->grav_field_f[p];
    uC[_MODEL_NVARS_*p+4] = (P*inv_gamma_m1)*(1.0/param->grav_field_g[p]) + (0.5*rho*(uvel*uvel+vvel*vvel+wvel*wvel))*param->grav_field_f[p];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  return(0);
}
