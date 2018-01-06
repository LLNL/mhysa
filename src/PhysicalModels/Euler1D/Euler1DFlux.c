/*! @file Euler1DFlux.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the hyperbolic flux for the 1D Euler equations over the domain.
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! 
  Compute the hyperbolic flux over the local domain.\n
*/
int Euler1DFlux(
                double  *f, /*!< Array to hold the computed flux (same size and layout as u) */
                double  *u, /*!< Array containing the conserved solution */
                int     dir,/*!< Spatial dimension (unused since this is a 1D system) */
                void    *s, /*!< Solver object of type #HyPar */
                double  t   /*!< Current time */
               )
{
  HyPar       *solver = (HyPar*)   s;
  Euler1D     *param  = (Euler1D*) solver->physics;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ns      = param->n_species;
  int nv      = param->n_vibeng;

  static int  index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double rho_s[ns], rho_t, v, E, E_v[nv], P, T;
    _Euler1DGetFlowVar_((u+nvars*p),rho_s,rho_t,v,E,E_v,P,T,param);
    _Euler1DSetFlux_((f+nvars*p),rho_s,rho_t,v,E,E_v,P,param);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
