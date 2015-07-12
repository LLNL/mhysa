/*! @file ShallowWater2DModifiedSolution.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute the modified solution for a balanced discretization scheme.
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the modified solution for the upwinding step in a balanced conservative
    finite-difference algorithm for the 2D shallow water equations. 
    \n\n
  Refer to:
  + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
    exact conservation property for the shallow water equations", Journal of
    Computational Physics, 208, 2005, pp. 206-227.
    http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater2DModifiedSolution(
                            double  *uC, /*!< The modified solution (same array size and layout as u) */
                            double  *u,  /*!< The solution (conserved variables) */
                            int     d,   /*!< Spatial dimension */
                            void    *s,  /*!< Solver object of type #HyPar */
                            void    *m,  /*!< MPI object of type #MPIVariables */
                            double  waqt /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*) s;
  ShallowWater2D  *param  = (ShallowWater2D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double h,uvel,vvel;
    _ShallowWater2DGetFlowVar_((u+_MODEL_NVARS_*p),h,uvel,vvel);
    uC[_MODEL_NVARS_*p+0] = h + param->b[p];
    uC[_MODEL_NVARS_*p+1] = h * uvel;
    uC[_MODEL_NVARS_*p+2] = h * vvel;
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  return(0);
}
