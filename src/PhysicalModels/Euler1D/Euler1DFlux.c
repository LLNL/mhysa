/*! @file Euler1DFlux.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the hyperbolic flux for the 1D Euler equations over the domain.
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! Compute the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ (e+p) u \end{array}\right]
    \f}
*/
int Euler1DFlux(
                double  *f, /*!< Array to hold the computed flux (same size and layout as u) */
                double  *u, /*!< Array containing the conserved solution */
                int     dir,/*!< Spatial dimension (unused since this is a 1D system) */
                void    *s, /*!< Solver object of type #HyPar */
                double  t   /*!< Current time */
               )
{
  HyPar             *solver = (HyPar*)   s;
  Euler1D           *param  = (Euler1D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  nvars   = _MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho, v, e, P;
    _Euler1DGetFlowVar_((u+nvars*p),rho,v,e,P,param);
    _Euler1DSetFlux_((f+nvars*p),rho,v,e,P);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

/*! Compute the stiff component of the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}_F\left({\bf u}\right) = A_f\left({\bf u}\right){\bf u}
    \f}
    where \f$A_f\left({\bf u}\right) \f$ is the fast Jacobian representing the 
    acoustic waves only. A linearized formulation is used where the fast Jacobian 
    \f$A_f\f$ is computed for the solution at the beginning of each time step in
    #Euler1DPreStep.
    \sa #_Euler1DSetStiffFlux_, #_Euler1DSetLinearizedStiffFlux_, #_Euler1DSetStiffJac_
*/
int Euler1DStiffFlux(
                double  *f, /*!< Array to hold the computed flux (same size and layout as u) */
                double  *u, /*!< Array containing the conserved solution */
                int     dir,/*!< Spatial dimension (unused since this is a 1D system) */
                void    *s, /*!< Solver object of type #HyPar */
                double  t   /*!< Current time */
               )
{
  HyPar             *solver = (HyPar*)   s;
  Euler1D           *param  = (Euler1D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  nvars   = _MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    _Euler1DSetLinearizedStiffFlux_((f+nvars*p),(u+nvars*p),(param->fast_jac+nvars*nvars*p));
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
