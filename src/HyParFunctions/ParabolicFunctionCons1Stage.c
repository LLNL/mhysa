/*! @file ParabolicFunctionCons1Stage.c
    @author Debojyoti Ghosh
    @brief Evaluate the parabolic term using a conservative spatial discretization
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Evaluate the parabolic term using a conservative finite-difference spatial discretization: 
    The parabolic term is assumed to be of the form:
    \f{equation}{
      {\bf P}\left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {\partial^2 {\bf g}_d\left(\bf u\right)} {\partial x_d^2},
    \f}
    which is discretized at a grid point as:
    \f{equation}{
      \left.{\bf P}\left({\bf u}\right)\right|_j = \sum_{d=0}^{D-1} \frac { \hat{\bf g}_{j+1/2} - \hat{\bf g}_{j-1/2} } {\Delta x_d^2},
    \f}
    where \f$d\f$ is the spatial dimension index, \f$D\f$ is the total number of spatial dimensions (#HyPar::ndims), and \f$j\f$ is
    the grid index along \f$d\f$. \f$\hat{\bf g}_d\f$ is the numerical approximation to the second primitive of 
    \f${\bf g}_d\left({\bf u}\right)\f$, computed using #HyPar::InterpolateInterfacesPar.

    \b Note: this form of the parabolic term \b does \b not allow for cross-derivatives.

    To use this form of the parabolic term:
    + specify \b "par_space_type" in solver.inp as \b "conservative-1stage" (#HyPar::spatial_type_par).
    + the physical model must specify \f${\bf g}_d\left({\bf u}\right)\f$ through #HyPar::GFunction.

    \b Reference: Liu, Y., Shu, C.-W., Zhang, M., "High order finite difference WENO schemes for nonlinear degenerate parabolic
                  equations", SIAM J. Sci. Comput., 33 (2), 2011, pp. 939-965, http://dx.doi.org/10.1137/100791002
*/
int ParabolicFunctionCons1Stage(
                                  double  *par, /*!< array to hold the computed parabolic term */
                                  double  *u,   /*!< solution */
                                  void    *s,   /*!< Solver object of type #HyPar */
                                  void    *m,   /*!< MPI object of type #MPIVariables */
                                  double  t     /*!< Current simulation time */
                               )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  double        *FluxI  = solver->fluxI; /* interface flux     array */
  double        *Func   = solver->fluxC; /* diffusion function array */
  int           d, v, i, done;
  _DECLARE_IERR_;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;
  int     size   = solver->npoints_local_wghosts;

  int index[ndims], index1[ndims], index2[ndims], dim_interface[ndims];

  _ArraySetValue_(par,size*nvars,0.0);
  if (!solver->GFunction) return(0); /* zero parabolic term */
  solver->count_par++;

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[d]++;
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];

    /* evaluate cell-centered flux */
    IERR solver->GFunction(Func,u,d,solver,t); CHECKERR(ierr);
    /* compute interface fluxes */
    IERR solver->InterpolateInterfacesPar(FluxI,Func,d,solver,mpi); CHECKERR(ierr);

    /* calculate the second derivative */
    done = 0; _ArraySetValue_(index,ndims,0);
    int p, p1, p2;
    while (!done) {
      _ArrayCopy1D_(index,index1,ndims);
      _ArrayCopy1D_(index,index2,ndims); index2[d]++;
      _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p);
      _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
      _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
      for (v=0; v<nvars; v++) 
        par[nvars*p+v] +=  ((dxinv[offset+ghosts+index[d]] * dxinv[offset+ghosts+index[d]]) 
                          * (FluxI[nvars*p2+v] - FluxI[nvars*p1+v]));
      _ArrayIncrementIndex_(ndims,dim,index,done);
    }

    offset += dim[d] + 2*ghosts;
  }

  if (solver->flag_ib) _ArrayBlockMultiply_(par,solver->iblank,size,nvars);
  return(0);
}
