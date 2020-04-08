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

/*! Roe's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R 
                         - \left| A\left({\bf u}_{j+1/2}^L,{\bf u}_{j+1/2}^R\right) \right|
                           \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    + Roe, P. L., “Approximate Riemann solvers, parameter vectors, and difference schemes,” Journal of
    Computational Physics, Vol. 43, No. 2, 1981, pp. 357–372, http://dx.doi.org/10.1016/0021-9991(81)90128-5.
    
    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when 
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces, 
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
int NavierStokes3DUpwindRoe(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension (x, y, or z) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*)    s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->physics;
  int             done;

  int *dim  = solver->dim_local;
  int nvars = solver->nvars;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], 
         L[nvars*nvars], DL[nvars*nvars], 
         modA[nvars*nvars];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->ghosts,pR);
      double udiff[nvars], uavg[nvars],udiss[nvars];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[nvars*p+0] - uL[nvars*p+0]);
      udiff[1] = 0.5 * (uR[nvars*p+1] - uL[nvars*p+1]);
      udiff[2] = 0.5 * (uR[nvars*p+2] - uL[nvars*p+2]);
      udiff[3] = 0.5 * (uR[nvars*p+3] - uL[nvars*p+3]);
      udiff[4] = 0.5 * (uR[nvars*p+4] - uL[nvars*p+4]);

      _NavierStokes3DRoeAverage_        (uavg,(u+nvars*pL),(u+nvars*pR),param);
      _NavierStokes3DEigenvalues_       (uavg,D,param,dir);
      _NavierStokes3DLeftEigenvectors_  (uavg,L,param,dir);
      _NavierStokes3DRightEigenvectors_ (uavg,R,param,dir);

      /* Harten's Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

      MatMult5(nvars,DL,D,L);
      MatMult5(nvars,modA,R,DL);
      MatVecMult5(nvars,udiss,modA,udiff);
      
      fI[nvars*p+0] = 0.5 * (fL[nvars*p+0]+fR[nvars*p+0]) - udiss[0];
      fI[nvars*p+1] = 0.5 * (fL[nvars*p+1]+fR[nvars*p+1]) - udiss[1];
      fI[nvars*p+2] = 0.5 * (fL[nvars*p+2]+fR[nvars*p+2]) - udiss[2];
      fI[nvars*p+3] = 0.5 * (fL[nvars*p+3]+fR[nvars*p+3]) - udiss[3];
      fI[nvars*p+4] = 0.5 * (fL[nvars*p+4]+fR[nvars*p+4]) - udiss[4];
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
int NavierStokes3DUpwindRF(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension (x, y, or z) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                          )
{
  HyPar           *solver = (HyPar*)    s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->physics;
  int             done,k;

  int *dim  = solver->dim_local;
  int nvars = solver->nvars;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[nvars], fcL[nvars], fcR[nvars], 
             ucL[nvars], ucR[nvars], fc[nvars];

      /* Roe-Fixed upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _NavierStokes3DLeftEigenvectors_(uavg,L,param,dir);
      _NavierStokes3DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(nvars,ucL,L,(uL+nvars*p));
      MatVecMult5(nvars,ucR,L,(uR+nvars*p));
      MatVecMult5(nvars,fcL,L,(fL+nvars*p));
      MatVecMult5(nvars,fcR,L,(fR+nvars*p));

      double eigL[nvars],eigC[nvars],eigR[nvars];
      _NavierStokes3DEigenvalues_((uL+nvars*p),D,param,dir); 
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _NavierStokes3DEigenvalues_((uR+nvars*p),D,param,dir); 
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _NavierStokes3DEigenvalues_(uavg,D,param,dir); 
      eigC[0] = D[0];
      eigC[1] = D[6];
      eigC[2] = D[12];
      eigC[3] = D[18];
      eigC[4] = D[24];

      for (k = 0; k < nvars; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult5(nvars,(fI+nvars*p),R,fc);
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

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when 
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces, 
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580
*/
int NavierStokes3DUpwindLLF(
                            double  *fI, /*!< Computed upwind interface flux */
                            double  *fL, /*!< Left-biased reconstructed interface flux */
                            double  *fR, /*!< Right-biased reconstructed interface flux */
                            double  *uL, /*!< Left-biased reconstructed interface solution */
                            double  *uR, /*!< Right-biased reconstructed interface solution */
                            double  *u,  /*!< Cell-centered solution */
                            int     dir, /*!< Spatial dimension (x, y, or z) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*)    s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->physics;
  int             done;

  int *dim  = solver->dim_local;
  int nvars = solver->nvars;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  double R[nvars*nvars], D[nvars*nvars], L[nvars*nvars];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[nvars], fcL[nvars], fcR[nvars], 
             ucL[nvars], ucR[nvars], fc[nvars];

      /* Local Lax-Friedrich upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,(uL+nvars*p),(uR+nvars*p),param);

      _NavierStokes3DLeftEigenvectors_(uavg,L,param,dir);
      _NavierStokes3DRightEigenvectors_(uavg,R,param,dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(nvars,ucL,L,(uL+nvars*p));
      MatVecMult5(nvars,ucR,L,(uR+nvars*p));
      MatVecMult5(nvars,fcL,L,(fL+nvars*p));
      MatVecMult5(nvars,fcR,L,(fR+nvars*p));

      double eigL[nvars],eigC[nvars],eigR[nvars];
      _NavierStokes3DEigenvalues_((uL+nvars*p),D,param,dir); 
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _NavierStokes3DEigenvalues_((uR+nvars*p),D,param,dir); 
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _NavierStokes3DEigenvalues_(uavg,D,param,dir); 
      eigC[0] = D[0];
      eigC[1] = D[6];
      eigC[2] = D[12];
      eigC[3] = D[18];
      eigC[4] = D[24];

      double alpha;
      alpha = max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));
      alpha = max3(absolute(eigL[4]),absolute(eigC[4]),absolute(eigR[4]));
      fc[4] = 0.5 * (fcL[4] + fcR[4] + alpha * (ucL[4]-ucR[4]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult5(nvars,(fI+nvars*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
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
