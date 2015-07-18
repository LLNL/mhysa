/*! @file ShallowWater1DTopography.c
    @author Debojyoti Ghosh
    @brief Contains the function to set the bottom topography
*/

#include <stdlib.h>
#include <math.h>
#include <arrayfunctions.h>
#include <io.h>
#include <physicalmodels/shallowwater1d.h>
#include <hypar.h>
#include <mpivars.h>

/*! Set the bottom topography over the domain - reads the topography
    data from a file, if available, else sets it to a constant
*/
int ShallowWater1DTopography(
                              void *s, /*!< Solver object of type #HyPar */
                              void *m  /*!< MPI object of type #MPIVariables */
                            )
{
  HyPar          *solver = (HyPar*) s;
  MPIVariables   *mpi    = (MPIVariables*) m;
  ShallowWater1D *param  = (ShallowWater1D*) solver->physics;
  double         *S      = param->b;
  int            d, done, *dim = solver->dim_local, 
                 ghosts = solver->ghosts;
  _DECLARE_IERR_;

  /* read topography from provided file, if available */
  IERR ReadArray(solver->ndims,1,solver->dim_global,solver->dim_local,solver->ghosts,
                 solver,mpi,NULL,S,"topography",&param->topo_flag); CHECKERR(ierr);
  if (!param->topo_flag) {
    /* if topography file not available, set it to zero */
    _ArraySetValue_(S,solver->npoints_local_wghosts,0.0);
  }

  /* if parallel, exchange MPI-boundary ghost point data */
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,1,solver->dim_local,
                                solver->ghosts,mpi,S); CHECKERR(ierr);


  if (param->bt_type) {
    /* if topography is periodic, then the overall problem must also be periodic
       (i.e. boundary conditions will be specified as periodic). Hence, 
       MPIExchangeBoundariesnD() will take care of setting the ghosts points 
       for multi-processor simulations. For single processor, set the ghost
       points accordingly. */
    int indexb[_MODEL_NDIMS_], indexi[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], 
        offset[_MODEL_NDIMS_];
    for (d = 0; d < _MODEL_NDIMS_; d++) {
      if (mpi->iproc[d] == 1) {
        _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = -ghosts;
        /* left boundary */
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = indexb[d] + dim[d] - ghosts;
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
        /* right boundary */
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_);
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
      }
    }
  } else {
    /* if topography is not periodic, extrapolate it at the boundaries */
    int indexb[_MODEL_NDIMS_], indexi[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], 
        offset[_MODEL_NDIMS_];
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
