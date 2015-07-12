/*! @file ShallowWater2DSourceUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the "upwind" source term at an interface (for a balanced finite-difference discretization of the 2D shallow water equations).
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <hypar.h>

/*! Compute the "upwind" source term in the balanced formulation introduced in the 
    reference below. The "upwind" state is just the arithmetic average of the left
    and right states.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater2DSourceUpwindLLF(
                                  double  *fI, /*!< Computed interface source term ("upwinded") */
                                  double  *fL, /*!< Left-biased interface source term */
                                  double  *fR, /*!< Right-biased interface source term */
                                  double  *u,  /*!< Solution (conserved variables) */
                                  int     dir, /*!< Spatial dimension */
                                  void    *s,  /*!< Solver object of type #HyPar */
                                  double  t    /*!< Current solution time */
                                 )
{
  HyPar     *solver = (HyPar*)    s;
  int       done,k;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (k = 0; k < _MODEL_NVARS_; k++) 
        (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Compute the "upwind" source term in the balanced formulation introduced in the 
    reference below. The "upwind" state is just the arithmetic average of the left
    and right states.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater2DSourceUpwindRoe(
                                  double  *fI, /*!< Computed interface source term ("upwinded") */
                                  double  *fL, /*!< Left-biased interface source term */
                                  double  *fR, /*!< Right-biased interface source term */
                                  double  *u,  /*!< Solution (conserved variables) */
                                  int     dir, /*!< Spatial dimension */
                                  void    *s,  /*!< Solver object of type #HyPar */
                                  double  t    /*!< Current solution time */
                                 )
{
  HyPar     *solver = (HyPar*)    s;
  int       done,k;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (k = 0; k < _MODEL_NVARS_; k++) 
        (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
