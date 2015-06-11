/*! @file Euler1DModifiedSolution.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute the modified solution for a balanced discretization scheme.
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the modified solution for the upwinding step in a balanced conservative
    finite-difference algorithm for the 1D Euler equations with gravitational sources. 
    If no gravitational forces exist, the modified solution is identical to the solution.
    \n\n
  Refer to 
  + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme 
    for the Gas Dynamics Equations Under Gravitational Fields", 
    Journal of Scientific Computing, 54, 2013, pp. 645-662,
    http://dx.doi.org/10.1007/s10915-012-9585-8. See Eq. (3.8) on 
    Page 651 on why this modification is needed.
*/
int Euler1DModifiedSolution(
                            double  *uC, /*!< The modified solution (same array size and layout as u) */
                            double  *u,  /*!< The solution (conserved variables) */
                            int     d,   /*!< Spatial dimension (unused since this is a 1D system) */
                            void    *s,  /*!< Solver object of type #HyPar */
                            void    *m,  /*!< MPI object of type #MPIVariables */
                            double  waqt /*!< Current solution time */
                           )
{
  HyPar         *solver = (HyPar*)         s;
  Euler1D       *param  = (Euler1D*)       solver->physics;

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
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p); p *= _MODEL_NVARS_;
    _ArrayScaleCopy1D_((u+p),(1.0/param->grav_field[p/_MODEL_NVARS_]),(uC+p),_MODEL_NVARS_);
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  return(0);
}
