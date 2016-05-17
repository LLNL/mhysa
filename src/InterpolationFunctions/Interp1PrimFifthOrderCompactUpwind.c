/*! @file Interp1PrimFifthOrderCompactUpwind.c
 *  @brief 5th order compact upwind scheme (component-wise application to vectors).
 *  @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <tridiagLU.h>
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

/*! @brief 5th order compact upwind reconstruction (component-wise) on a uniform grid

    Computes the interpolated values of the first primitive of a function \f${\bf f}\left({\bf u}\right)\f$
    at the interfaces from the cell-centered values of the function using the fifth order compact upwind scheme on a 
    uniform grid. The first primitive is defined as a function \f${\bf h}\left({\bf u}\right)\f$ that satisfies:
    \f{equation}{
      {\bf f}\left({\bf u}\left(x\right)\right) = \frac{1}{\Delta x} \int_{x-\Delta x/2}^{x+\Delta x/2} {\bf h}\left({\bf u}\left(\zeta\right)\right)d\zeta,
    \f}
    where \f$x\f$ is the spatial coordinate along the dimension of the interpolation. This function computes the 5th order compact upwind numerical approximation 
    \f$\hat{\bf f}_{j+1/2} \approx {\bf h}_{j+1/2}\f$ as:
    \f{align}{
      \frac{3}{10}\hat{\bf f}_{j-1/2} + \frac{6}{10}\hat{\bf f}_{j+1/2} + \frac{1}{10}\hat{\bf f}_{j+3/2} = \frac{1}{30}{\bf f}_{j-1} + \frac{19}{30}{\bf f}_j + \frac{1}{3}{\bf f}_{j+1}.
    \f}
    The resulting tridiagonal system is solved using tridiagLU() (see also #TridiagLU, tridiagLU.h).

    \b Implementation \b Notes:
    + This method assumes a uniform grid in the spatial dimension corresponding to the interpolation.
    + The method described above corresponds to a left-biased interpolation. The corresponding right-biased
      interpolation can be obtained by reflecting the equations about interface j+1/2.
    + The scalar interpolation method is applied to the vector function in a component-wise manner.
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
    + Ghosh, D., Baeder, J. D., Compact Reconstruction Schemes with Weighted ENO Limiting 
      for Hyperbolic Conservation Laws, SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706,
      http://dx.doi.org/10.1137/110857659
 */
int Interp1PrimFifthOrderCompactUpwind(
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
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  CompactScheme   *compact= (CompactScheme*)  solver->compact;
  TridiagLU       *lu     = (TridiagLU*)      solver->lusolver;
  int             sys,Nsys,d,v;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int *stride= solver->stride_with_ghosts;

  /* define some constants */
  static const double one_third            = 1.0/3.0,
                      thirteen_by_sixty    = 13.0/60.0,
                      fortyseven_by_sixty  = 47.0/60.0,
                      twentyseven_by_sixty = 27.0/60.0,
                      one_by_twenty        = 1.0/20.0,
                      one_by_thirty        = 1.0/30.0,
                      nineteen_by_thirty   = 19.0/30.0,
                      three_by_ten         = 3.0/10.0,
                      six_by_ten           = 6.0/10.0,
                      one_by_ten           = 1.0/10.0;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  /* calculate total number of tridiagonal systems to solve */
  _ArrayProduct1D_(bounds_outer,ndims,Nsys); Nsys *= nvars;

  /* Allocate arrays for tridiagonal system */
  double *A = compact->A;
  double *B = compact->B;
  double *C = compact->C;
  double *R = compact->R;

#pragma omp parallel for schedule(auto) default(shared) private(sys,d,index_outer,indexC,indexI)
  for (sys=0; sys < N_outer; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      if (upw > 0) {
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        qm3 = qm1 - 2*stride[dir];
        qm2 = qm1 -   stride[dir];
        qp1 = qm1 +   stride[dir];
        qp2 = qm1 + 2*stride[dir];
      } else {
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        qm3 = qm1 + 2*stride[dir];
        qm2 = qm1 +   stride[dir];
        qp1 = qm1 -   stride[dir];
        qp2 = qm1 - 2*stride[dir];
      }

      /* Defining stencil points */
      double *fm3, *fm2, *fm1, *fp1, *fp2;
      fm3 = fC+qm3*nvars;
      fm2 = fC+qm2*nvars;
      fm1 = fC+qm1*nvars;
      fp1 = fC+qp1*nvars;
      fp2 = fC+qp2*nvars;

      if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
          || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {

        /* Use 5th order upwind at the physical boundaries */
        _ArraySetValue_((A+Nsys*indexI[dir]+sys*nvars),nvars,0.0)
        _ArraySetValue_((B+Nsys*indexI[dir]+sys*nvars),nvars,1.0)
        _ArraySetValue_((C+Nsys*indexI[dir]+sys*nvars),nvars,0.0)
        for (v=0; v<nvars; v++) {
          (R+Nsys*indexI[dir]+sys*nvars)[v] =   one_by_thirty         * fm3[v]
                                              - thirteen_by_sixty     * fm2[v]
                                              + fortyseven_by_sixty   * fm1[v]
                                              + twentyseven_by_sixty  * fp1[v]
                                              - one_by_twenty         * fp2[v];
        }

      } else {

        /* Use 5th order upwind at the physical boundaries */
        if (upw > 0) {
          _ArraySetValue_((A+Nsys*indexI[dir]+sys*nvars),nvars,three_by_ten);
          _ArraySetValue_((B+Nsys*indexI[dir]+sys*nvars),nvars,six_by_ten  );
          _ArraySetValue_((C+Nsys*indexI[dir]+sys*nvars),nvars,one_by_ten  );
        } else {
          _ArraySetValue_((C+Nsys*indexI[dir]+sys*nvars),nvars,three_by_ten);
          _ArraySetValue_((B+Nsys*indexI[dir]+sys*nvars),nvars,six_by_ten  );
          _ArraySetValue_((A+Nsys*indexI[dir]+sys*nvars),nvars,one_by_ten  );
        }
        for (v=0; v<nvars; v++) {
          (R+Nsys*indexI[dir]+sys*nvars)[v] =   one_by_thirty      * fm2[v]
                                              + nineteen_by_thirty * fm1[v]
                                              + one_third          * fp1[v];
        }
      }
    }
  }

#ifdef serial

  /* Solve the tridiagonal system */
  IERR tridiagLU(A,B,C,R,dim[dir]+1,Nsys,lu,NULL); CHECKERR(ierr);

#else

  /* Solve the tridiagonal system */
  /* all processes except the last will solve without the last interface to avoid overlap */
  if (mpi->ip[dir] != mpi->iproc[dir]-1)  { IERR tridiagLU(A,B,C,R,dim[dir]  ,Nsys,lu,&mpi->comm[dir]); CHECKERR(ierr); }
  else                                    { IERR tridiagLU(A,B,C,R,dim[dir]+1,Nsys,lu,&mpi->comm[dir]); CHECKERR(ierr); }

  /* Now get the solution to the last interface from the next proc */
  double *sendbuf = compact->sendbuf;
  double *recvbuf = compact->recvbuf;
  MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
  if (mpi->ip[dir]) for (d=0; d<Nsys; d++) sendbuf[d] = R[d];
  if (mpi->ip[dir] != mpi->iproc[dir]-1) MPI_Irecv(recvbuf,Nsys,MPI_DOUBLE,mpi->ip[dir]+1,214,mpi->comm[dir],&req[0]);
  if (mpi->ip[dir])                      MPI_Isend(sendbuf,Nsys,MPI_DOUBLE,mpi->ip[dir]-1,214,mpi->comm[dir],&req[1]);
  MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
  if (mpi->ip[dir] != mpi->iproc[dir]-1) for (d=0; d<Nsys; d++) R[d+Nsys*dim[dir]] = recvbuf[d];

#endif

  /* save the solution to fI */
#pragma omp parallel for schedule(auto) default(shared) private(sys,d,index_outer,indexC,indexI)
  for (sys=0; sys < N_outer; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      _ArrayCopy1D_((R+sys*nvars+Nsys*indexI[dir]),(fI+nvars*p),nvars);
    }
  }

  return(0);
}
