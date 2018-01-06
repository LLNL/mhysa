/*! @file NavierStokes3DUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the upwind flux at grid interfaces for the 3D Navier Stokes equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*! Rusanov's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R 
                         - \max_{j,j+1} \nu_j \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    where \f$\nu = c + \left|u\right|\f$.
    + Rusanov, V. V., "The calculation of the interaction of non-stationary shock waves and obstacles," USSR 
    Computational Mathematics and Mathematical Physics, Vol. 1, No. 2, 1962, pp. 304–320
*/
int NavierStokes3DUpwindRusanov(
                                double  *fI, /*!< Computed upwind interface flux */
                                double  *fL, /*!< Left-biased reconstructed interface flux */
                                double  *fR, /*!< Right-biased reconstructed interface flux */
                                double  *uL, /*!< Left-biased reconstructed interface solution */
                                double  *uR, /*!< Right-biased reconstructed interface solution */
                                double  *u,  /*!< Cell-centered solution */
                                int     dir, /*!< Spatial dimension (x,y, or z) */
                                void    *s,  /*!< Solver object of type #HyPar */
                                double  t    /*!< Current solution time */
                               )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             done, k;
  _DECLARE_IERR_;

  int nvars = solver->nvars;
  int ns    = param->n_species;
  int nv    = param->n_vibeng;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  static int  bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_],
              index_outer[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_], 
              indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_];

  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[2] = dim[2]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[dir]++;

  double udiff[nvars],uavg[nvars];

  done = 0; _ArraySetValue_(index_outer,_MODEL_NDIMS_,0);
  while (!done) {

    _ArrayCopy1D_(index_outer,index_inter,_MODEL_NDIMS_);

    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {

      int p; _ArrayIndex1D_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,ghosts,pR);

      /* Modified Rusanov's upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,(u+nvars*pL),(u+nvars*pR),param);
      for (k = 0; k < nvars; k++) udiff[k] = 0.5 * (uR[nvars*p+k] - uL[nvars*p+k]);

      double rho_s[ns], rho_t, vel[_MODEL_NDIMS_], E, E_v[nv], P, T, c;

      _NavierStokes3DGetFlowVar_( (u+nvars*pL),rho_s,rho_t,
                                  vel[_XDIR_],vel[_YDIR_],vel[_ZDIR_],
                                  E,E_v,P,T,param );
      _NavierStokes3DSpeedOfSound_(c,(param->gamma),P,rho_t);
      double alphaL = c + absolute(vel[dir]);

      _NavierStokes3DGetFlowVar_( (u+nvars*pR),rho_s,rho_t,
                                  vel[_XDIR_],vel[_YDIR_],vel[_ZDIR_],
                                  E,E_v,P,T,param );
      _NavierStokes3DSpeedOfSound_(c,(param->gamma),P,rho_t);
      double alphaR = c + absolute(vel[dir]);

      _NavierStokes3DGetFlowVar_( uavg,rho_s,rho_t,
                                  vel[_XDIR_],vel[_YDIR_],vel[_ZDIR_],
                                  E,E_v,P,T,param );
      _NavierStokes3DSpeedOfSound_(c,(param->gamma),P,rho_t);
      double alphaavg = c + absolute(vel[dir]);

      double alpha  = max3(alphaL,alphaR,alphaavg);

      for (k = 0; k < nvars; k++) {
        fI[nvars*p+k] = 0.5 * (fL[nvars*p+k]+fR[nvars*p+k]) - alpha*udiff[k];
      }

    }

    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);

  }

  return(0);
}
