/*! @file Euler1DComputeCFL.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute maximum CFL over the domain for the 1D Euler equations physical model.
*/

#include <stdlib.h>
#include <basic.h>
#include <math.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! Computes the maximum CFL number over the domain. Note that the CFL
    is computed over the local domain on this processor only.
*/
double Euler1DComputeCFL(
                          void    *s, /*!< Solver object of type #HyPar */
                          void    *m, /*!< MPI object of type #MPIVariables */
                          double  dt, /*!< Time step size for which to compute the CFL */
                          double  t   /*!< Time */
                        )
{
  HyPar     *solver = (HyPar*)   s;
  Euler1D   *param  = (Euler1D*) solver->physics;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;

  int ns      = param->n_species;
  int nv      = param->n_vibeng;

  int index[ndims];
  double *u   = solver->u;

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {

    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double rho_s[ns], rho_t, v, E, E_v[nv], P, T, c, dxinv, local_cfl;

    _Euler1DGetFlowVar_((u+nvars*p),rho_s,rho_t,v,E,E_v,P,T,param);
    _Euler1DSpeedOfSound_(c,(param->gamma),P,rho_t);
    _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */

    local_cfl = (absolute(v)+c)*dt*dxinv; /* local cfl for this grid point */
    if (local_cfl > max_cfl) max_cfl = local_cfl;

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
