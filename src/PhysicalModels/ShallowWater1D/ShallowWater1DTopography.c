/*! @file ShallowWater1DTopography.c
    @author Debojyoti Ghosh
    @brief Contains the function to set the bottom topography
*/

#include <stdlib.h>
#include <math.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <hypar.h>
#include <mpivars.h>

/*! Set the bottom topography over the domain. */
int ShallowWater1DTopography(
                              void *s, /*!< Solver object of type #HyPar */
                              void *m  /*!< MPI object of type #MPIVariables */
                            )
{
  HyPar          *solver = (HyPar*) s;
  MPIVariables   *mpi    = (MPIVariables*) m;
  ShallowWater1D *param  = (ShallowWater1D*) solver->physics;

  double  *S      = param->b;
  int     *dim    = solver->dim_local;
  int     ghosts  = solver->ghosts;
  int     index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_],
          offset[_MODEL_NDIMS_], d, done;

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_);
  for (d=0; d<_MODEL_NDIMS_; d++) bounds[d] += 2*ghosts;
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  /* set the value of the topography */
  done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p;         _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double xcoord; _GetCoordinate_(_XDIR_,index[_XDIR_]-ghosts,dim,ghosts,solver->x,xcoord);
    if (param->bt_type == 0) {
      S[p] = 0.0;
    } else if (param->bt_type == 1) {
      S[p] = 5.0 * exp(-0.4*(xcoord-5.0)*(xcoord-5.0));
    } else if (param->bt_type == 2) {
      if (xcoord < 4.0)      S[p] = 0.0;
      else if (xcoord < 8.0) S[p] = 4.0;
      else                   S[p] = 0.0;
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  if (param->bt_type != -1) {
    /* if topography is not periodic, extrapolate it at the boundaries */
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
