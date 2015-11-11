/*! @file Interp1PrimFifthOrderHCWENO.c
    @author Debojyoti Ghosh
    @brief hybrid compact-WENO5 Scheme (Component-wise application to vectors)
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

/*! @brief 5th order hybrid compact-WENO reconstruction (component-wise) on a uniform grid

    Computes the interpolated values of the first primitive of a function \f${\bf f}\left({\bf u}\right)\f$
    at the interfaces from the cell-centered values of the function using the fifth order hybrid compact-WENO scheme on a 
    uniform grid. The tridiagonal system is solved using tridiagLU() (see also #TridiagLU, tridiagLU.h). See references 
    below for a complete description of the method implemented here.

    \b Implementation \b Notes:
    + This method assumes a uniform grid in the spatial dimension corresponding to the interpolation.
    + The scalar interpolation method is applied to the vector function in a component-wise manner.
    + The WENO weights are computed in WENOFifthOrderCalculateWeights().
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
    + Pirozzoli, S., Conservative Hybrid Compact-WENO Schemes for Shock-Turbulence Interaction, J. Comput. Phys., 178 (1), 2002, pp. 81-117, http://dx.doi.org/10.1006/jcph.2002.7021
    + Ren, Y.-X., Liu, M., Zhang, H., A characteristic-wise hybrid compact-WENO scheme for solving hyperbolic conservation laws, J. Comput. Phys., 192 (2), 2003, pp. 365-386, 
      http://dx.doi.org/10.1016/j.jcp.2003.07.006
 */
int Interp1PrimFifthOrderHCWENO(
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
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  TridiagLU       *lu     = (TridiagLU*)      solver->lusolver;
  int             sys,Nsys,d;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  static const double one_half           = 1.0/2.0;
  static const double one_third          = 1.0/3.0;
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

  /* calculate total number of tridiagonal systems to solve */
  _ArrayProduct1D_(bounds_outer,ndims,Nsys); Nsys *= nvars;

  /* Allocate arrays for tridiagonal system */
  double *A = weno->A;
  double *B = weno->B;
  double *C = weno->C;
  double *R = weno->R;

#pragma omp parallel for schedule(auto) default(shared) private(sys,d,index_outer,indexC,indexI)
  for (sys=0; sys < N_outer; sys++) {
    _ArrayIndexnD_(ndims,sys,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
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
      int v; 
      for (v=0; v<nvars; v++)  {
        /* Defining stencil points */
        double fm3, fm2, fm1, fp1, fp2;
        fm3 = fC[qm3*nvars+v];
        fm2 = fC[qm2*nvars+v];
        fm1 = fC[qm1*nvars+v];
        fp1 = fC[qp1*nvars+v];
        fp2 = fC[qp2*nvars+v];

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

        /* calculate the hybridization parameter */
        double sigma;
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Standard WENO at physical boundaries */
          sigma = 0.0;
        } else {
          double cuckoo, df_jm12, df_jp12, df_jp32, r_j, r_jp1, r_int;
          cuckoo = (0.9*weno->rc / (1.0-0.9*weno->rc)) * weno->xi * weno->xi;
          df_jm12 = fm1 - fm2;
          df_jp12 = fp1 - fm1;
          df_jp32 = fp2 - fp1;
          r_j   = (absolute(2*df_jp12*df_jm12)+cuckoo)/(df_jp12*df_jp12+df_jm12*df_jm12+cuckoo);
          r_jp1 = (absolute(2*df_jp32*df_jp12)+cuckoo)/(df_jp32*df_jp32+df_jp12*df_jp12+cuckoo);
          r_int = min(r_j, r_jp1);
          sigma = min((r_int/weno->rc), 1.0); 
        }

        if (upw > 0) {
          A[sys*nvars+v+Nsys*indexI[dir]] = one_half  * sigma;
          B[sys*nvars+v+Nsys*indexI[dir]] = 1.0;
          C[sys*nvars+v+Nsys*indexI[dir]] = one_sixth * sigma;
        } else {
          C[sys*nvars+v+Nsys*indexI[dir]] = one_half  * sigma;
          B[sys*nvars+v+Nsys*indexI[dir]] = 1.0;
          A[sys*nvars+v+Nsys*indexI[dir]] = one_sixth * sigma;
        }

        double fWENO, fCompact;
        fWENO    = w1*f1 + w2*f2 + w3*f3;
        fCompact = one_sixth * (one_third*fm2 + 19.0*one_third*fm1 + 10.0*one_third*fp1);
        R[sys*nvars+v+Nsys*indexI[dir]] =   sigma*fCompact + (1.0-sigma)*fWENO;
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
  double *sendbuf = weno->sendbuf;
  double *recvbuf = weno->recvbuf;
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
