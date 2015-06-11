/*! @file Euler1DGravityField.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute the gravitational field
*/

#include <stdlib.h>
#include <math.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>
#include <mpivars.h>

/*! Compute the gravitational field over the domain. See
    + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
      for the Gas Dynamics Equations Under Gravitational Fields",
      Journal of Scientific Computing, 54, 2013, pp. 645-662.
      http://dx.doi.org/10.1007/s10915-012-9585-8\n
    for details on the treatment of the gravitational source term.
    This function computes the graviational potential function. If
    there is no gravity, the gravitational field is 1.0 all over 
    the domain.
*/
int Euler1DGravityField(
                        void *s, /*!< Solver object of type #HyPar */
                        void *m  /*!< MPI object of type #MPIVariables */
                       )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  Euler1D       *param  = (Euler1D*)      solver->physics;

  double  *S      = param->grav_field;
  int     *dim    = solver->dim_local;
  int     ghosts  = solver->ghosts;
  int     index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_],
          offset[_MODEL_NDIMS_], d, done;

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_);
  for (d=0; d<_MODEL_NDIMS_; d++) bounds[d] += 2*ghosts;
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  /* set the value of the gravity field */
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p;                _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double xcoord;        _GetCoordinate_(_XDIR_,index[_XDIR_]-ghosts,dim,ghosts,solver->x,xcoord);
    if (param->grav_type == 0) {
      S[p] = exp(-param->grav*xcoord);
    } else if (param->grav_type == 1) {
      double pi = 4.0 * atan(1.0);
      double phi = -sin(2*pi*xcoord)/(2*pi);
      S[p] = exp(-phi);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  if (param->grav_type != 1) {
    /* a sensible simulation will not specify peridic boundary conditions 
    * along a direction in which gravity acts (unless the gravitational field
    * itself is sinusoidal), so extrapolate the gravity field at the boundaries */
    int indexb[_MODEL_NDIMS_], indexi[_MODEL_NDIMS_];
    for (d = 0; d < _MODEL_NDIMS_; d++) {
      /* left boundary */
      if (!mpi->ip[d]) {
        _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = -ghosts;
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = ghosts-1-indexb[d];
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
      }
      /* right boundary */
      if (mpi->ip[d] == mpi->iproc[d]-1) {
        _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = dim[d];
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = dim[d]-1-indexb[d];
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
      }
    }
  }

  return(0);
}
