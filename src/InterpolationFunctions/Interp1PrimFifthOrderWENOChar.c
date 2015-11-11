/*! @file Interp1PrimFifthOrderWENOChar.c
    @author Debojyoti Ghosh
    @brief Characteristic-based WENO5 Scheme
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_omp
#include <omp.h>
#endif

#undef  _MINIMUM_GHOSTS_
/*! \def _MINIMUM_GHOSTS_
 * Minimum number of ghost points required for this interpolation 
 * method.
*/
#define _MINIMUM_GHOSTS_ 3

/*! @brief 5th order WENO reconstruction (characteristic-based) on a uniform grid

    Computes the interpolated values of the first primitive of a function \f${\bf f}\left({\bf u}\right)\f$
    at the interfaces from the cell-centered values of the function using the fifth order WENO scheme on a 
    uniform grid. The first primitive is defined as a function \f${\bf h}\left({\bf u}\right)\f$ that satisfies:
    \f{equation}{
      {\bf f}\left({\bf u}\left(x\right)\right) = \frac{1}{\Delta x} \int_{x-\Delta x/2}^{x+\Delta x/2} {\bf h}\left({\bf u}\left(\zeta\right)\right)d\zeta,
    \f}
    where \f$x\f$ is the spatial coordinate along the dimension of the interpolation. This function computes the 5th order WENO numerical approximation \f$\hat{\bf f}_{j+1/2} \approx {\bf h}_{j+1/2}\f$ as the convex combination of three 3rd order methods:
    \f{align}{
        &\ \omega_1\ \times\ \left[ \hat{\alpha}^k_{j+1/2} = \frac{1}{3} {\alpha}^k_{j-2} - \frac{7}{6} {\alpha}^k_{j-1} + \frac{11}{6} {\alpha}^k_j \right]\\
      + &\ \omega_2\ \times\ \left[ \hat{\alpha}^k_{j+1/2} = -\frac{1}{6} {\alpha}^k_{j-1} + \frac{5}{6} {\alpha}^k_j + \frac{1}{3} {\alpha}^k_{j+1} \right]\\
      + &\ \omega_3\ \times\ \left[ \hat{\alpha}^k_{j+1/2} = \frac{1}{3} {\alpha}^k_j + \frac{5}{6} {\alpha}^k_{j+1} - \frac{1}{6} {\alpha}^k_{j+2}  \right]\\
      \Rightarrow &\ \hat{\alpha}^k_{j+1/2} = \frac{\omega_1}{3} {\alpha}^k_{j-2} - \frac{1}{6}(7\omega_1+\omega_2){\alpha}^k_{j-1} + \frac{1}{6}(11\omega_1+5\omega_2+2\omega_3){\alpha}^k_j + \frac{1}{6}(2\omega_2+5\omega_3){\alpha}^k_{j+1} - \frac{\omega_3}{6}{\alpha}^k_{j+2},
    \f}
    where
    \f{equation}{
      \alpha^k = {\bf l}_k \cdot {\bf f},\ k=1,\cdots,n
    \f}
    is the \f$k\f$-th characteristic quantity, and \f${\bf l}_k\f$ is the \f$k\f$-th left eigenvector, \f${\bf r}_k\f$ is the \f$k\f$-th right eigenvector, and \f$n\f$ is #HyPar::nvars. The nonlinear weights \f$\omega_k; k=1,2,3\f$ are the WENO weights computed in WENOFifthOrderCalculateWeightsChar(). The final interpolated function is computed from the interpolated characteristic quantities as:
    \f{equation}{
      \hat{\bf f}_{j+1/2} = \sum_{k=1}^n \alpha^k_{j+1/2} {\bf r}_k
    \f}

    \b Implementation \b Notes:
    + This method assumes a uniform grid in the spatial dimension corresponding to the interpolation.
    + The method described above corresponds to a left-biased interpolation. The corresponding right-biased
      interpolation can be obtained by reflecting the equations about interface j+1/2.
    + The left and right eigenvectors are computed at an averaged quantity at j+1/2. Thus, this function requires
      functions to compute the average state, and the left and right eigenvectors. These are provided by the physical
      model through
      - #HyPar::GetLeftEigenvectors() 
      - #HyPar::GetRightEigenvectors()
      - #HyPar::AveragingFunction() 

      If these functions are not provided by the physical model, then a characteristic-based interpolation cannot be used.
    + The function computes the interpolant for the entire grid in one call. It loops over all the grid lines along the interpolation direction
      and carries out the 1D interpolation along these grid lines.
    + Location of cell-centers and cell interfaces along the spatial dimension of the interpolation is shown in the following figure:
      @image html chap1_1Ddomain.png
      @image latex chap1_1Ddomain.eps width=0.9\textwidth

    \b Function \b arguments:

    Argument  | Type      | Explanation             
    --------- | --------- | ---------------------------------------------
    fI        | double*   | Array to hold the computed interpolant at the grid interfaces. This array must have the same layout as the solution, but with \b no \b ghost \b points. Its size should be the same as u in all dimensions, except dir (the dimension along which to interpolate) along which it should be larger by 1 (number of interfaces is 1 more than the number of interior cell centers).
    fC        | double*   | Array with the cell-centered values of the flux function \f${\bf f}\left({\bf u}\right)\f$. This array must have the same layout and size as the solution, \b with \b ghost \b points. 
    u         | double*   | The solution array \f${\bf u}\f$ (with ghost points). If the interpolation is characteristic based, this is needed to compute the eigendecomposition. For a multidimensional problem, the layout is as follows: u is a contiguous 1D array of size (nvars*dim[0]*dim[1]*...*dim[D-1]) corresponding to the multi-dimensional solution, with the following ordering - nvars, dim[0], dim[1], ..., dim[D-1], where nvars is the number of solution components (#HyPar::nvars), dim is the local size (#HyPar::dim_local), D is the number of spatial dimensions.
    x         | double*   | The grid array (with ghost points). This is used only by non-uniform-grid interpolation methods. For multidimensional problems, the layout is as follows: x is a contiguous 1D array of size (dim[0]+dim[1]+...+dim[D-1]), with the spatial coordinates along dim[0] stored from 0,...,dim[0]-1, the spatial coordinates along dim[1] stored along dim[0],...,dim[0]+dim[1]-1, and so forth.
    upw       | int       | Upwinding direction: if positive, a left-biased interpolant will be computed; if negative, a right-biased interpolant will be computed. If the interpolation method is central, then this has no effect.
    dir       | int       | Spatial dimension along which to interpolate (eg: 0 for 1D; 0 or 1 for 2D; 0,1 or 2 for 3D)
    s         | void*     | Solver object of type #HyPar: the following variables are needed - #HyPar::ghosts, #HyPar::ndims, #HyPar::nvars, #HyPar::dim_local.
    m         | void*     | MPI object of type #MPIVariables: this is needed only by compact interpolation method that need to solve a global implicit system across MPI ranks.
    uflag     | int       | A flag indicating if the function being interpolated \f${\bf f}\f$ is the solution itself \f${\bf u}\f$ (if 1, \f${\bf f}\left({\bf u}\right) \equiv {\bf u}\f$).


    \b Reference: 
    + Jiang, G.-S., Shu, C.-W., Efficient Implementation of Weighted ENO Schemes, J. Comput. Phys., 126 (1), 1996, pp. 202-228, http://dx.doi.org/10.1006/jcph.1996.0130
 */
int Interp1PrimFifthOrderWENOChar(
                                    double *fI,  /*!< Array of interpolated function values at the interfaces */
                                    double *fC,  /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
                                    double *u,   /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
                                    double *x,   /*!< Grid coordinates */
                                    int    upw,  /*!< Upwind direction (left or right biased) */
                                    int    dir,  /*!< Spatial dimension along which to interpolation */
                                    void   *s,   /*!< Object of type #HyPar containing solver-related variables */
                                    void   *m,   /*!< Object of type #MPIVariables containing MPI-related variables */
                                    int    uflag /*!< Flag to indicate if \f$f(u) \equiv u\f$, i.e, if the solution is being reconstructed */
                                 )
{
  HyPar           *solver = (HyPar*) s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  int             i, k, v;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;

  double *ww1, *ww2, *ww3;
  ww1 = weno->w1 + (upw < 0 ? 2*weno->size : 0) + (uflag ? weno->size : 0) + weno->offset[dir];
  ww2 = weno->w2 + (upw < 0 ? 2*weno->size : 0) + (uflag ? weno->size : 0) + weno->offset[dir];
  ww3 = weno->w3 + (upw < 0 ? 2*weno->size : 0) + (uflag ? weno->size : 0) + weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double R[nvars*nvars], L[nvars*nvars], uavg[nvars], fchar[nvars];

#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      /* 1D indices of the stencil grid points */
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }

      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&u[nvars*qm1],&u[nvars*qp1],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double fm3, fm2, fm1, fp1, fp2;
        fm3 = fm2 = fm1 = fp1 = fp2 = 0;
        for (k = 0; k < nvars; k++) {
          fm3 += L[v*nvars+k] * fC[qm3*nvars+k];
          fm2 += L[v*nvars+k] * fC[qm2*nvars+k];
          fm1 += L[v*nvars+k] * fC[qm1*nvars+k];
          fp1 += L[v*nvars+k] * fC[qp1*nvars+k];
          fp2 += L[v*nvars+k] * fC[qp2*nvars+k];
        }

        /* Candidate stencils and their optimal weights*/
        double f1, f2, f3;
        f1 = (2*one_sixth)*fm3 - (7.0*one_sixth)*fm2 + (11.0*one_sixth)*fm1;
        f2 = (-one_sixth)*fm2 + (5.0*one_sixth)*fm1 + (2*one_sixth)*fp1;
        f3 = (2*one_sixth)*fm1 + (5*one_sixth)*fp1 - (one_sixth)*fp2;

        /* calculate WENO weights */
        double w1,w2,w3;
        w1 = *(ww1+p*nvars+v);
        w2 = *(ww2+p*nvars+v);
        w3 = *(ww3+p*nvars+v);

        /* fifth order WENO approximation of the characteristic flux */
        fchar[v] = w1*f1 + w2*f2 + w3*f3;

      }

      /* calculate the interface u from the characteristic u */
      IERR MatVecMult(nvars,(fI+nvars*p),R,fchar); CHECKERR(ierr);

    }
  }

  return(0);
}
