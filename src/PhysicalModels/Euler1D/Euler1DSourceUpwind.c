/*! @file Euler1DSourceUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the "upwind" source term at an interface (for a balanced finite-difference discretization of the 1D Euler equations with gravitational source terms).
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! Compute the "upwind" source term in the balanced formulation introduced in the 
    reference below. The "upwind" state is just the arithmetic average of the left
    and right states.
    + Xing, Shu, "High Order Well-Balanced WENO Scheme for the Gas Dynamics Equations 
                  Under Gravitational Fields", J. Sci. Comput., 54, 2013, pp. 645--662,
                  http://dx.doi.org/10.1007/s10915-012-9585-8.
*/
int Euler1DSourceUpwindLLF(
                            double  *fI, /*!< Computed interface source term ("upwinded") */
                            double  *fL, /*!< Left-biased interface source term */
                            double  *fR, /*!< Right-biased interface source term */
                            double  *u,  /*!< Solution (conserved variables) */
                            int     dir, /*!< Spatial dimension (unused since this is a 1D case) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                          )
{
  HyPar     *solver = (HyPar*)    s;
  int       done,k;
  _DECLARE_IERR_;

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
      /* Local Lax-Friedrich upwinding scheme */
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
    + Xing, Shu, "High Order Well-Balanced WENO Scheme for the Gas Dynamics Equations 
                  Under Gravitational Fields", J. Sci. Comput., 54, 2013, pp. 645--662,
                  http://dx.doi.org/10.1007/s10915-012-9585-8.
*/
int Euler1DSourceUpwindRoe(
                            double  *fI, /*!< Computed interface source term ("upwinded") */
                            double  *fL, /*!< Left-biased interface source term */
                            double  *fR, /*!< Right-biased interface source term */
                            double  *u,  /*!< Solution (conserved variables) */
                            int     dir, /*!< Spatial dimension (unused since this is a 1D case) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            double  t    /*!< Current solution time */
                          )
{
  HyPar     *solver = (HyPar*)    s;
  int       done,k;
  _DECLARE_IERR_;

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
      /* Local Lax-Friedrich upwinding scheme */
      for (k = 0; k < _MODEL_NVARS_; k++) 
        (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
