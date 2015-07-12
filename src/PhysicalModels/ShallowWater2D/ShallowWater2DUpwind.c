/*! @file ShallowWater2DUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the upwind flux at grid interfaces for the 2D shallow water equations.
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <hypar.h>

/*! Roe's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R 
                         - \left| A\left({\bf u}_{j+1/2}^L,{\bf u}_{j+1/2}^R\right) \right|
                           \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    + Roe, P. L., “Approximate Riemann solvers, parameter vectors, and difference schemes,” Journal of
    Computational Physics, Vol. 43, No. 2, 1981, pp. 357–372, http://dx.doi.org/10.1016/0021-9991(81)90128-5.
    
    This upwinding scheme is modified for the balanced discretization of the 2D shallow water equations. See the following reference:
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006

*/
int ShallowWater2DUpwindRoe(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*) s;
  ShallowWater2D  *param  = (ShallowWater2D*) solver->physics;
  int             done,k;

  int ndims = solver->ndims;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_], 
                DL[_MODEL_NVARS_*_MODEL_NVARS_], modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);

      _ShallowWater2DRoeAverage_         (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param); 
      _ShallowWater2DEigenvalues_        (uavg,D,param,dir);
      _ShallowWater2DLeftEigenvectors_   (uavg,L,param,dir);
      _ShallowWater2DRightEigenvectors_  (uavg,R,param,dir);

      k = 0; D[k] = absolute(D[k]);
      k = 4; D[k] = absolute(D[k]);
      k = 8; D[k] = absolute(D[k]);

      MatMult3(_MODEL_NVARS_,DL,D,L);
      MatMult3(_MODEL_NVARS_,modA,R,DL);

      udiss[0] = modA[0*_MODEL_NVARS_+0]*udiff[0] + modA[0*_MODEL_NVARS_+1]*udiff[1] + modA[0*_MODEL_NVARS_+2]*udiff[2];
      udiss[1] = modA[1*_MODEL_NVARS_+0]*udiff[0] + modA[1*_MODEL_NVARS_+1]*udiff[1] + modA[1*_MODEL_NVARS_+2]*udiff[2];
      udiss[2] = modA[2*_MODEL_NVARS_+0]*udiff[0] + modA[2*_MODEL_NVARS_+1]*udiff[1] + modA[2*_MODEL_NVARS_+2]*udiff[2];

      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Characteristic-based local Lax-Friedrich upwinding scheme.
    \f{align}{
      \alpha_{j+1/2}^{k,L} &= \sum_{k=1}^2 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,L}, \\
      \alpha_{j+1/2}^{k,R} &= \sum_{k=1}^2 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,R}, \\
      v_{j+1/2}^{k,L} &= \sum_{k=1}^2 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,L}, \\
      v_{j+1/2}^{k,R} &= \sum_{k=1}^2 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,R}, \\
      \alpha_{j+1/2}^k &= \frac{1}{2}\left[ \alpha_{j+1/2}^{k,L} + \alpha_{j+1/2}^{k,R} - \left(\max_{\left[j,j+1\right]} \lambda\right) \left( v_{j+1/2}^{k,R} - v_{j+1/2}^{k,L} \right) \right], \\
      {\bf f}_{j+1/2} &= \sum_{k=1}^2 \alpha_{j+1/2}^k {\bf r}_{j+1/2}^k
    \f}
    where \f${\bf l}\f$, \f${\bf r}\f$, and \f$\lambda\f$ are the left-eigenvectors, right-eigenvectors and eigenvalues. The subscripts denote the grid locations.
    + C.-W. Shu, and S. Osher, "Efficient implementation of essentially non-oscillatory schemes, II", J. Comput. Phys., 83 (1989), pp. 32–78, http://dx.doi.org/10.1016/0021-9991(89)90222-2.

    This upwinding scheme is modified for the balanced discretization of the 1D shallow water equations. See the following reference:
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater2DUpwindLLF(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*) s;
  ShallowWater2D  *param  = (ShallowWater2D*) solver->physics;
  int             done,k;

  int ndims   = solver->ndims;
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Local Lax-Friedrich upwinding scheme */

      _ShallowWater2DRoeAverage_       (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _ShallowWater2DEigenvalues_      (uavg,D,param,dir);
      _ShallowWater2DLeftEigenvectors_ (uavg,L,param,dir);
      _ShallowWater2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult3(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _ShallowWater2DEigenvalues_((u+_MODEL_NVARS_*pL),D,param,dir);
        eigL = D[k*_MODEL_NVARS_+k];
        _ShallowWater2DEigenvalues_((u+_MODEL_NVARS_*pR),D,param,dir);
        eigR = D[k*_MODEL_NVARS_+k];
        _ShallowWater2DEigenvalues_(uavg,D,param,dir);
        eigC = D[k*_MODEL_NVARS_+k];

        double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult3(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

