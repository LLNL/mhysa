/*! @file NavierStokes2DUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the upwind flux at grid interfaces for the 2D Navier Stokes equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mathfunctions.h>
#include <hypar.h>

/*! Roe's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R 
                         - \left| A\left({\bf u}_{j+1/2}^L,{\bf u}_{j+1/2}^R\right) \right|
                           \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    + Roe, P. L., “Approximate Riemann solvers, parameter vectors, and difference schemes,” Journal of
    Computational Physics, Vol. 43, No. 2, 1981, pp. 357–372, http://dx.doi.org/10.1016/0021-9991(81)90128-5.
    
    This upwinding scheme is modified for the balanced discretization of the 2D Navier Stokes equations when 
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces, 
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted

*/
int NavierStokes2DUpwindRoe(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension (x or y) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int *dim  = solver->dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      _NavierStokes2DRoeAverage_        (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

       /* Harten's Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      double kappa = max(param->grav_field_g[pL],param->grav_field_g[pR]);
      k=0;  D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=5;  D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=10; D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=15; D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss,modA,udiff);
      
      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5 * (fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Characteristic-based Roe-fixed upwinding scheme.
    \f{align}{
      \alpha_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,L}, \\
      \alpha_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,R}, \\
      v_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,L}, \\
      v_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,R}, \\
      \alpha_{j+1/2}^k &= \left\{ \begin{array}{cc} \alpha_{j+1/2}^{k,L} & {\rm if}\ \lambda_{j,j+1/2,j+1} > 0 \\ \alpha_{j+1/2}^{k,R} & {\rm if}\ \lambda_{j,j+1/2,j+1} < 0 \\ \frac{1}{2}\left[ \alpha_{j+1/2}^{k,L} + \alpha_{j+1/2}^{k,R} - \left(\max_{\left[j,j+1\right]} \lambda\right) \left( v_{j+1/2}^{k,R} - v_{j+1/2}^{k,L} \right) \right] & {\rm otherwise} \end{array}\right., \\
      {\bf f}_{j+1/2} &= \sum_{k=1}^3 \alpha_{j+1/2}^k {\bf r}_{j+1/2}^k
    \f}
    where \f${\bf l}\f$, \f${\bf r}\f$, and \f$\lambda\f$ are the left-eigenvectors, right-eigenvectors and eigenvalues. The subscripts denote the grid locations.
    + C.-W. Shu, and S. Osher, "Efficient implementation of essentially non-oscillatory schemes, II", J. Comput. Phys., 83 (1989), pp. 32–78, http://dx.doi.org/10.1016/0021-9991(89)90222-2.

    Note that this upwinding scheme cannot be used for solving flows with non-zero gravitational forces.
*/
int NavierStokes2DUpwindRF(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension (x or y) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar    *solver = (HyPar*)    s;
  NavierStokes2D  *param  = (NavierStokes2D*)  solver->physics;
  int      done,k;

  int *dim  = solver->dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _NavierStokes2DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_ (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _NavierStokes2DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,dir);
      eigL[0] = D[0];
      eigL[1] = D[5];
      eigL[2] = D[10];
      eigL[3] = D[15];
      _NavierStokes2DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,dir);
      eigR[0] = D[0];
      eigR[1] = D[5];
      eigR[2] = D[10];
      eigR[3] = D[15];
      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      eigC[0] = D[0];
      eigC[1] = D[5];
      eigC[2] = D[10];
      eigC[3] = D[15];

      for (k = 0; k < _MODEL_NVARS_; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Characteristic-based local Lax-Friedrich upwinding scheme.
    \f{align}{
      \alpha_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,L}, \\
      \alpha_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,R}, \\
      v_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,L}, \\
      v_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,R}, \\
      \alpha_{j+1/2}^k &= \frac{1}{2}\left[ \alpha_{j+1/2}^{k,L} + \alpha_{j+1/2}^{k,R} - \left(\max_{\left[j,j+1\right]} \lambda\right) \left( v_{j+1/2}^{k,R} - v_{j+1/2}^{k,L} \right) \right], \\
      {\bf f}_{j+1/2} &= \sum_{k=1}^3 \alpha_{j+1/2}^k {\bf r}_{j+1/2}^k
    \f}
    where \f${\bf l}\f$, \f${\bf r}\f$, and \f$\lambda\f$ are the left-eigenvectors, right-eigenvectors and eigenvalues. The subscripts denote the grid locations.
    + C.-W. Shu, and S. Osher, "Efficient implementation of essentially non-oscillatory schemes, II", J. Comput. Phys., 83 (1989), pp. 32–78, http://dx.doi.org/10.1016/0021-9991(89)90222-2.

    This upwinding scheme is modified for the balanced discretization of the 2D Navier Stokes equations when 
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces, 
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted
*/
int NavierStokes2DUpwindLLF(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension (x or y) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar    *solver = (HyPar*)    s;
  NavierStokes2D  *param  = (NavierStokes2D*)  solver->physics;
  int      done;

  int *dim  = solver->dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];
      double kappa = max(param->grav_field_g[pL],param->grav_field_g[pR]);

      /* Local Lax-Friedrich upwinding scheme */

      _NavierStokes2DRoeAverage_        (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _NavierStokes2DEigenvalues_((u+_MODEL_NVARS_*pL),D,param,dir);
      eigL[0] = D[0];
      eigL[1] = D[5];
      eigL[2] = D[10];
      eigL[3] = D[15];
      _NavierStokes2DEigenvalues_((u+_MODEL_NVARS_*pR),D,param,dir);
      eigR[0] = D[0];
      eigR[1] = D[5];
      eigR[2] = D[10];
      eigR[3] = D[15];
      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      eigC[0] = D[0];
      eigC[1] = D[5];
      eigC[2] = D[10];
      eigC[3] = D[15];

      double alpha;
      alpha = kappa * max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Steger-Warming Flux-Splitting scheme
    + Steger, J.L., Warming, R.F., "Flux vector splitting of the inviscid gasdynamic equations with 
      application to finite-difference methods", J. Comput. Phys., 40(2), 1981, pp. 263-293,
      http://dx.doi.org/10.1016/0021-9991(81)90210-2.

    Note that this method cannot be used for flows with non-zero gravitational forces.
*/
int NavierStokes2DUpwindSWFS(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                            )
{
  HyPar             *solver = (HyPar*)    s;
  NavierStokes2D    *param  = (NavierStokes2D*)  solver->physics;
  int               done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double fp[_MODEL_NVARS_], fm[_MODEL_NVARS_],uavg[_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double rho,vx,vy,e,P,c,gamma=param->gamma,term,Mach,lp[_MODEL_NVARS_],lm[_MODEL_NVARS_];

      /* Steger Warming flux splitting */
      _NavierStokes2DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param); 
      _NavierStokes2DGetFlowVar_(uavg,rho,vx,vy,e,P,param);
      Mach = (dir==_XDIR_ ? vx : vy) / sqrt(gamma*P/rho);

      if (Mach < -1.0) {

        _ArrayCopy1D_((fR+_MODEL_NVARS_*p),(fI+_MODEL_NVARS_*p),_MODEL_NVARS_);

      } else if (Mach < 1.0) {

        double kx = 0, ky = 0;
        kx = (dir==_XDIR_ ? 1.0 : 0.0);
        ky = (dir==_YDIR_ ? 1.0 : 0.0);

        _NavierStokes2DGetFlowVar_((uL+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);
        c = sqrt(gamma*P/rho);
        term = rho/(2.0*gamma);
        lp[0] = lp[1] = kx*vx + ky*vy;
        lp[2] = lp[0] + c;
        lp[3] = lp[0] - c;
        for (k=0; k<_MODEL_NVARS_; k++) if (lp[k] < 0.0) lp[k] = 0.0;

        fp[0] = term * (2.0*(gamma-1.0)*lp[0] + lp[2] + lp[3]);
        fp[1] = term * (2.0*(gamma-1.0)*lp[0]*vx + lp[2]*(vx+c*kx) + lp[3]*(vx-c*kx));
        fp[2] = term * (2.0*(gamma-1.0)*lp[0]*vy + lp[2]*(vy+c*ky) + lp[3]*(vy-c*ky));
        fp[3] = term * ((gamma-1.0)*lp[0]*(vx*vx+vy*vy) + 0.5*lp[2]*((vx+c*kx)*(vx+c*kx) + (vy+c*ky)*(vy+c*ky)) 
                        + 0.5*lp[3]*((vx-c*kx)*(vx-c*kx) + (vy-c*ky)*(vy-c*ky)) 
                        + ((3.0-gamma)*(lp[2]+lp[3])*c*c)/(2.0*(gamma-1.0)) );

        _NavierStokes2DGetFlowVar_((uR+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);
        c = sqrt(gamma*P/rho);
        term = rho/(2.0*gamma);
        lm[0] = lm[1] = kx*vx + ky*vy;
        lm[2] = lm[0] + c;
        lm[3] = lm[0] - c;
        for (k=0; k<_MODEL_NVARS_; k++) if (lm[k] > 0.0) lm[k] = 0.0;

        fm[0] = term * (2.0*(gamma-1.0)*lm[0] + lm[2] + lm[3]);
        fm[1] = term * (2.0*(gamma-1.0)*lm[0]*vx + lm[2]*(vx+c*kx) + lm[3]*(vx-c*kx));
        fm[2] = term * (2.0*(gamma-1.0)*lm[0]*vy + lm[2]*(vy+c*ky) + lm[3]*(vy-c*ky));
        fm[3] = term * ((gamma-1.0)*lm[0]*(vx*vx+vy*vy) + 0.5*lm[2]*((vx+c*kx)*(vx+c*kx) + (vy+c*ky)*(vy+c*ky)) 
                        + 0.5*lm[3]*((vx-c*kx)*(vx-c*kx) + (vy-c*ky)*(vy-c*ky)) 
                        + ((3.0-gamma)*(lm[2]+lm[3])*c*c)/(2.0*(gamma-1.0)) );

        _ArrayAdd1D_((fI+_MODEL_NVARS_*p),fp,fm,_MODEL_NVARS_);

      } else {

        _ArrayCopy1D_((fL+_MODEL_NVARS_*p),(fI+_MODEL_NVARS_*p),_MODEL_NVARS_);

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Rusanov's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R 
                         - \max_{j,j+1} \nu_j \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    where \f$\nu = c + \left|u\right|\f$.
    + Rusanov, V. V., "The calculation of the interaction of non-stationary shock waves and obstacles," USSR 
    Computational Mathematics and Mathematical Physics, Vol. 1, No. 2, 1962, pp. 304–320
    
    This upwinding scheme is modified for the balanced discretization of the 2D Navier Stokes equations when 
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces, 
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted

*/
int NavierStokes2DUpwindRusanov(
                                double  *fI, /*!< Computed upwind interface flux */
                                double  *fL, /*!< Left-biased reconstructed interface flux */
                                double  *fR, /*!< Right-biased reconstructed interface flux */
                                double  *uL, /*!< Left-biased reconstructed interface solution */
                                double  *uR, /*!< Right-biased reconstructed interface solution */
                                double  *u,  /*!< Cell-centered solution */
                                int     dir, /*!< Spatial dimension (x or y) */
                                void    *s,  /*!< Solver object of type #HyPar */
                                double  t    /*!< Current solution time */
                               )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int *dim  = solver->dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_];

      /* Modified Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      _NavierStokes2DRoeAverage_        (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

      double c, vel[_MODEL_NDIMS_], rho,E,P;
      _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pL),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaL = c + absolute(vel[dir]), betaL = absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pR),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaR = c + absolute(vel[dir]), betaR = absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_(uavg,rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaavg = c + absolute(vel[dir]), betaavg = absolute(vel[dir]);

      double kappa  = max(param->grav_field_g[pL],param->grav_field_g[pR]);
      double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

      fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Modified Rusanov's upwinding scheme: NavierStokes2DUpwindRusanov() modified as described in the
    following paper (for consistent characteristic-based splitting):
    + Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of Atmospheric Flows with 
      Characteristic-Based Flux Partitioning", Submitted (http://arxiv.org/abs/1510.05751).

*/
int NavierStokes2DUpwindRusanovModified(
                                          double  *fI, /*!< Computed upwind interface flux */
                                          double  *fL, /*!< Left-biased reconstructed interface flux */
                                          double  *fR, /*!< Right-biased reconstructed interface flux */
                                          double  *uL, /*!< Left-biased reconstructed interface solution */
                                          double  *uR, /*!< Right-biased reconstructed interface solution */
                                          double  *u,  /*!< Cell-centered solution */
                                          int     dir, /*!< Spatial dimension (x or y) */
                                          void    *s,  /*!< Solver object of type #HyPar */
                                          double  t    /*!< Current solution time */
                                        )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int *dim  = solver->dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Modified Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      _NavierStokes2DRoeAverage_        (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

      double c, vel[_MODEL_NDIMS_], rho,E,P;
      _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pL),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaL = c + absolute(vel[dir]), betaL = absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pR),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaR = c + absolute(vel[dir]), betaR = absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_(uavg,rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaavg = c + absolute(vel[dir]), betaavg = absolute(vel[dir]);

      double kappa  = max(param->grav_field_g[pL],param->grav_field_g[pR]);
      double alpha  = kappa*max3(alphaL,alphaR,alphaavg);
      double beta   = kappa*max3(betaL,betaR,betaavg);

      _ArraySetValue_(D,_MODEL_NVARS_,0.0);
      D[0]  = alpha;
      D[5]  = (dir == _XDIR_ ? alpha : beta);
      D[10] = (dir == _YDIR_ ? alpha : beta);
      D[15] = beta;
      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss,modA,udiff);

      fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0])-udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1])-udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2])-udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3])-udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The Roe upwinding scheme (#NavierStokes2DUpwindRoe) for the partitioned hyperbolic flux that comprises
    of the acoustic waves only (see #NavierStokes2DStiffFlux, #_NavierStokes2DSetStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$ u\pm a\f$ are used.
*/
int NavierStokes2DUpwinddFRoe(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                             )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      _NavierStokes2DRoeAverage_        (uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

       /* Harten's Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      double kappa = max(param->grav_field_g[pL],param->grav_field_g[pR]);
      k=0;  D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=5;  D[k] = (dir == _YDIR_ ? 0.0 : kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) ) );
      k=10; D[k] = (dir == _XDIR_ ? 0.0 : kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) ) );
      k=15; D[k] = 0.0;

      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss,modA,udiff);
      
      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5 * (fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The characteristic-based Roe-fixed upwinding scheme (#NavierStokes2DUpwindRF) for the partitioned hyperbolic flux 
    that comprises of the acoustic waves only (see #NavierStokes2DStiffFlux, #_NavierStokes2DSetStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$ u\pm a\f$ are used.
*/
int NavierStokes2DUpwinddFRF(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                            )
{
  HyPar    *solver = (HyPar*)    s;
  NavierStokes2D  *param  = (NavierStokes2D*)  solver->physics;
  int      done,k;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _NavierStokes2DRoeAverage_(uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);

      _NavierStokes2DEigenvalues_      (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_ (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pL),D,param,dir);
      eigL[0] = D[0];
      eigL[1] = (dir == _YDIR_ ? 0.0 : D[5]);
      eigL[2] = (dir == _XDIR_ ? 0.0 : D[10]);
      eigL[3] = 0.0;
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pR),D,param,dir);
      eigR[0] = D[0];
      eigR[1] = (dir == _YDIR_ ? 0.0 : D[5]);
      eigR[2] = (dir == _XDIR_ ? 0.0 : D[10]);
      eigR[3] = 0.0;
      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      eigC[0] = D[0];
      eigC[1] = (dir == _YDIR_ ? 0.0 : D[5]);
      eigC[2] = (dir == _XDIR_ ? 0.0 : D[10]);
      eigC[3] = 0.0;

      for (k = 0; k < _MODEL_NVARS_; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The characteristic-based local Lax-Friedrich upwinding scheme (#NavierStokes2DUpwindLLF) for the partitioned hyperbolic 
    flux that comprises of the acoustic waves only (see #NavierStokes2DStiffFlux, #_NavierStokes2DSetStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$ u\pm a\f$ are used.
*/
int NavierStokes2DUpwinddFLLF(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                             )
{
  HyPar    *solver = (HyPar*)    s;
  NavierStokes2D  *param  = (NavierStokes2D*)  solver->physics;
  int      done;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];
      double kappa = max(param->grav_field_g[pL],param->grav_field_g[pR]);

      /* Local Lax-Friedrich upwinding scheme */

      _NavierStokes2DRoeAverage_        (uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pL),D,param,dir);
      eigL[0] = D[0];
      eigL[1] = (dir == _YDIR_ ? 0.0 : D[5]);
      eigL[2] = (dir == _XDIR_ ? 0.0 : D[10]);
      eigL[3] = 0.0;
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pR),D,param,dir);
      eigR[0] = D[0];
      eigR[1] = (dir == _YDIR_ ? 0.0 : D[5]);
      eigR[2] = (dir == _XDIR_ ? 0.0 : D[10]);
      eigR[3] = 0.0;
      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      eigC[0] = D[0];
      eigC[1] = (dir == _YDIR_ ? 0.0 : D[5]);
      eigC[2] = (dir == _XDIR_ ? 0.0 : D[10]);
      eigC[3] = 0.0;

      double alpha;
      alpha = kappa * max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The modified Rusanov upwinding scheme (NavierStokes2DUpwindRusanovModified()) for the partitioned hyperbolic flux that comprises
    of the acoustic waves only (see #NavierStokes2DStiffFlux, #_NavierStokes2DSetStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$ u\pm a\f$ are used.
*/
int NavierStokes2DUpwinddFRusanovModified(
                                            double  *fI, /*!< Computed upwind interface flux */
                                            double  *fL, /*!< Left-biased reconstructed interface flux */
                                            double  *fR, /*!< Right-biased reconstructed interface flux */
                                            double  *uL, /*!< Left-biased reconstructed interface solution */
                                            double  *uR, /*!< Right-biased reconstructed interface solution */
                                            double  *u,  /*!< Cell-centered solution */
                                            int     dir, /*!< Spatial dimension (x or y) */
                                            void    *s,  /*!< Solver object of type #HyPar */
                                            double  t    /*!< Current solution time */
                                         )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[_MODEL_NVARS_], udiss[_MODEL_NVARS_], uavg[_MODEL_NVARS_];

      /* Rusanov's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      _NavierStokes2DRoeAverage_        (uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

      double c, vel[_MODEL_NDIMS_], rho,E,P;
      _NavierStokes2DGetFlowVar_((uref+_MODEL_NVARS_*pL),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaL = c + absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_((uref+_MODEL_NVARS_*pR),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaR = c + absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_(uavg,rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      double alphaavg = c + absolute(vel[dir]);

      double kappa  = max(param->grav_field_g[pL],param->grav_field_g[pR]);
      double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

      _ArraySetValue_(D,_MODEL_NVARS_*_MODEL_NVARS_,0.0);
      D[0]  = alpha;
      D[5]  = (dir == _YDIR_ ? 0.0 : alpha);
      D[10] = (dir == _XDIR_ ? 0.0 : alpha);
      D[15] = 0.0;

      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss,modA,udiff);

      fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The Roe upwinding scheme (#NavierStokes2DUpwindRoe) for the partitioned hyperbolic flux that comprises
    of the entropy waves only (see #NavierStokes2DNonStiffFlux, #_NavierStokes2DSetStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$u\f$ are used.
*/
int NavierStokes2DUpwindFdFRoe(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                             )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      int k;
      double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_],
             udiss_total[_MODEL_NVARS_],udiss_stiff[_MODEL_NVARS_];
      double delta = 0.000001, delta2 = delta*delta;
      double kappa = max(param->grav_field_g[pL],param->grav_field_g[pR]);

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      /* Compute total dissipation */
      _NavierStokes2DRoeAverage_        (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);
      k=0;  D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=5;  D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=10; D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=15; D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss_total,modA,udiff);

      /* Compute dissipation corresponding to acoustic modes */
      _NavierStokes2DRoeAverage_        (uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);
      k=0;  D[k] = kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=5;  D[k] = (dir == _YDIR_ ? 0.0 : kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) ) );
      k=10; D[k] = (dir == _XDIR_ ? 0.0 : kappa * (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) ) );
      k=15; D[k] = 0.0;
      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss_stiff,modA,udiff);
     
     /* Compute the dissipation term for the entropy modes */
      _ArraySubtract1D_(udiss,udiss_total,udiss_stiff,_MODEL_NVARS_);

      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5 * (fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The characteristic-based Roe-fixed upwinding scheme (#NavierStokes2DUpwindRF) for the partitioned hyperbolic flux 
    that comprises of the entropy waves only (see #NavierStokes2DNonStiffFlux, #_NavierStokes2DSetNonStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$u\f$ are used.
*/
int NavierStokes2DUpwindFdFRF(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                            )
{
  HyPar    *solver = (HyPar*)    s;
  NavierStokes2D  *param  = (NavierStokes2D*)  solver->physics;
  int      done,k;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _NavierStokes2DRoeAverage_(uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);

      _NavierStokes2DEigenvalues_      (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_ (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pL),D,param,dir);
      eigL[0] = 0.0;
      eigL[1] = (dir == _XDIR_ ? 0.0 : D[5]);
      eigL[2] = (dir == _YDIR_ ? 0.0 : D[10]);
      eigL[3] = D[15];
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pR),D,param,dir);
      eigR[0] = 0.0; 
      eigR[1] = (dir == _XDIR_ ? 0.0 : D[5]);
      eigR[2] = (dir == _YDIR_ ? 0.0 : D[10]);
      eigR[3] = D[15];
      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      eigC[0] = 0.0; 
      eigC[1] = (dir == _XDIR_ ? 0.0 : D[5]);
      eigC[2] = (dir == _YDIR_ ? 0.0 : D[10]);
      eigC[3] = D[15];

      for (k = 0; k < _MODEL_NVARS_; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The characteristic-based local Lax-Friedrich upwinding scheme (#NavierStokes2DUpwindLLF) for the partitioned hyperbolic 
    flux that comprises of the entropy waves only (see #NavierStokes2DNonStiffFlux, #_NavierStokes2DSetNonStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$u\f$ are used.
*/
int NavierStokes2DUpwindFdFLLF(
                              double  *fI, /*!< Computed upwind interface flux */
                              double  *fL, /*!< Left-biased reconstructed interface flux */
                              double  *fR, /*!< Right-biased reconstructed interface flux */
                              double  *uL, /*!< Left-biased reconstructed interface solution */
                              double  *uR, /*!< Right-biased reconstructed interface solution */
                              double  *u,  /*!< Cell-centered solution */
                              int     dir, /*!< Spatial dimension (x or y) */
                              void    *s,  /*!< Solver object of type #HyPar */
                              double  t    /*!< Current solution time */
                             )
{
  HyPar    *solver = (HyPar*)    s;
  NavierStokes2D  *param  = (NavierStokes2D*)  solver->physics;
  int      done;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];
      double kappa = max(param->grav_field_g[pL],param->grav_field_g[pR]);

      /* Local Lax-Friedrich upwinding scheme */

      _NavierStokes2DRoeAverage_        (uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);
      _NavierStokes2DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pL),D,param,dir);
      eigL[0] = 0.0;
      eigL[1] = (dir == _XDIR_ ? 0.0 : D[5]);
      eigL[2] = (dir == _YDIR_ ? 0.0 : D[10]);
      eigL[3] = D[15];
      _NavierStokes2DEigenvalues_((uref+_MODEL_NVARS_*pR),D,param,dir);
      eigR[0] = 0.0; 
      eigR[1] = (dir == _XDIR_ ? 0.0 : D[5]);
      eigR[2] = (dir == _YDIR_ ? 0.0 : D[10]);
      eigR[3] = D[15];
      _NavierStokes2DEigenvalues_(uavg,D,param,dir);
      eigC[0] = 0.0; 
      eigC[1] = (dir == _XDIR_ ? 0.0 : D[5]);
      eigC[2] = (dir == _YDIR_ ? 0.0 : D[10]);
      eigC[3] = D[15];

      double alpha;
      alpha = kappa * max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The modified Rusanov upwinding scheme (NavierStokes2DUpwindRusanovModified()) for the partitioned hyperbolic flux that comprises
    of the entropy waves only (see #NavierStokes2DNonStiffFlux, #_NavierStokes2DSetNonStiffFlux_). Thus, only the 
    characteristic fields / eigen-modes corresponding to \f$u\f$ are used.
*/
int NavierStokes2DUpwindFdFRusanovModified(
                                            double  *fI, /*!< Computed upwind interface flux */
                                            double  *fL, /*!< Left-biased reconstructed interface flux */
                                            double  *fR, /*!< Right-biased reconstructed interface flux */
                                            double  *uL, /*!< Left-biased reconstructed interface solution */
                                            double  *uR, /*!< Right-biased reconstructed interface solution */
                                            double  *u,  /*!< Cell-centered solution */
                                            int     dir, /*!< Spatial dimension (x or y) */
                                            void    *s,  /*!< Solver object of type #HyPar */
                                            double  t    /*!< Current solution time */
                                          )
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             done;

  int     *dim  = solver->dim_local;
  double  *uref = param->solution;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_], 
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[_MODEL_NVARS_], udiss[_MODEL_NVARS_], uavg[_MODEL_NVARS_],
             udiss_total[_MODEL_NVARS_],udiss_acoustic[_MODEL_NVARS_];

      /* Modified Rusanov's upwinding scheme */
      double c, vel[_MODEL_NDIMS_], rho,E,P, alphaL, alphaR, alphaavg,
             betaL, betaR, betaavg, alpha, beta,
             kappa  = max(param->grav_field_g[pL],param->grav_field_g[pR]);


      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

      /* Compute total dissipation */
      _NavierStokes2DRoeAverage_        (uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),param);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);
      _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pL),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      alphaL = c + absolute(vel[dir]); 
      betaL = absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pR),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      alphaR = c + absolute(vel[dir]); 
      betaR = absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_(uavg,rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      alphaavg = c + absolute(vel[dir]); 
      betaavg = absolute(vel[dir]);
      alpha  = kappa*max3(alphaL,alphaR,alphaavg);
      beta   = kappa*max3(betaL,betaR,betaavg);
      _ArraySetValue_(D,_MODEL_NVARS_,0.0);
      D[0]  = alpha;
      D[5]  = (dir == _XDIR_ ? alpha : beta);
      D[10] = (dir == _YDIR_ ? alpha : beta);
      D[15] = beta;
      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss_total,modA,udiff);

      /* Compute dissipation for the linearized acoustic modes */
      _NavierStokes2DRoeAverage_        (uavg,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param);
      _NavierStokes2DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes2DRightEigenvectors_ (uavg,R,param,dir);
      _NavierStokes2DGetFlowVar_((uref+_MODEL_NVARS_*pL),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      alphaL = c + absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_((uref+_MODEL_NVARS_*pR),rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      alphaR = c + absolute(vel[dir]);
      _NavierStokes2DGetFlowVar_(uavg,rho,vel[0],vel[1],E,P,param);
      c = sqrt(param->gamma*P/rho);
      alphaavg = c + absolute(vel[dir]);
      kappa  = max(param->grav_field_g[pL],param->grav_field_g[pR]);
      alpha  = kappa*max3(alphaL,alphaR,alphaavg);
      _ArraySetValue_(D,_MODEL_NVARS_*_MODEL_NVARS_,0.0);
      D[0]  = alpha;
      D[5]  = (dir == _YDIR_ ? 0.0 : alpha);
      D[10] = (dir == _XDIR_ ? 0.0 : alpha);
      D[15] = 0.0;
      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss_acoustic,modA,udiff);

      /* Compute dissipation for the entropy modes */
      _ArraySubtract1D_(udiss,udiss_total,udiss_acoustic,_MODEL_NVARS_);

      fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
      fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3]) - udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}
