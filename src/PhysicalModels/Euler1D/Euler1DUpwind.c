/*! @file Euler1DUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the upwind flux at grid interfaces for the 1D Euler equations.
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
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
int Euler1DUpwindRusanov(
                          double  *fI, /*!< Computed upwind interface flux */
                          double  *fL, /*!< Left-biased reconstructed interface flux */
                          double  *fR, /*!< Right-biased reconstructed interface flux */
                          double  *uL, /*!< Left-biased reconstructed interface solution */
                          double  *uR, /*!< Right-biased reconstructed interface solution */
                          double  *u,  /*!< Cell-centered solution */
                          int     dir, /*!< Spatial dimension (unused since this is a 1D system) */
                          void    *s,  /*!< Solver object of type #HyPar */
                          double  t    /*!< Current solution time */
                        )
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int ns    = param->n_species;
  int nv    = param->n_vibeng;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  double udiff[nvars], uavg[nvars];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {

    _ArrayCopy1D_(index_outer,index_inter,ndims);

    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {

      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

      _Euler1DRoeAverage_(uavg,(u+nvars*pL),(u+nvars*pR),param);
      for (k = 0; k < nvars; k++) udiff[k] = 0.5 * (uR[nvars*p+k] - uL[nvars*p+k]);

      double rho_s[ns], rho_t, uvel, E, E_v[nv], P, T, c;

      _Euler1DGetFlowVar_((u+nvars*pL),rho_s,rho_t,uvel,E,E_v,P,T,param);
      _Euler1DSpeedOfSound_(c,(param->gamma),P,rho_t);
      double alphaL = c + absolute(uvel);

      _Euler1DGetFlowVar_((u+nvars*pR),rho_s,rho_t,uvel,E,E_v,P,T,param);
      _Euler1DSpeedOfSound_(c,(param->gamma),P,rho_t);
      double alphaR = c + absolute(uvel);

      _Euler1DGetFlowVar_(uavg,rho_s,rho_t,uvel,E,E_v,P,T,param);
      _Euler1DSpeedOfSound_(c,(param->gamma),P,rho_t);
      double alphaavg = c + absolute(uvel);

      double alpha  = max3(alphaL,alphaR,alphaavg);

      for (k = 0; k < nvars; k++) {
        fI[nvars*p+k] = 0.5 * (fL[nvars*p+k]+fR[nvars*p+k]) - alpha*udiff[k];
      }

    }

    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
