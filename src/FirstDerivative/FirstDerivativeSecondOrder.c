/*! @file FirstDerivativeSecondOrder.c
    @author Debojyoti Ghosh
    @brief Second order finite-difference approximation to first derivative
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <firstderivative.h>

#include <mpivars.h>
#include <hypar.h>
typedef MPIVariables  MPIContext;
typedef HyPar         SolverContext;

#ifdef with_omp
#include <omp.h>
#endif

/*! Computes the second-order finite-difference approximation to the first derivative 
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} \frac{1}{2}\left(-3f_i+4f_{i+1}-f_{i+2}\right) & i = -g \\ \frac{1}{2}\left( f_{i+1} - f_{i-1} \right) & -g+1 \leq i \leq N+g-2 \\ \frac{1}{2}\left( f_{i-2} -4f_{i-1}+3f_i \right) & i = N+g-1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative, \f$g\f$ is the number of ghost points, and \f$N\f$ is the number of grid points (not including the ghost points) in the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + The first derivative is computed at the ghost points too. Thus, biased schemes are used
      at and near the boundaries.
    + \b Df and \b f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.
*/
int FirstDerivativeSecondOrderCentral(
                                        double  *Df,  /*!< Array to hold the computed first derivative (with ghost points) */
                                        double  *f,   /*!< Array containing the grid point function values whose first 
                                                           derivative is to be computed (with ghost points) */
                                        int     dir,  /*!< The spatial dimension along which the derivative is computed */
                                        int     bias, /*!< Forward or backward differencing for non-central 
                                                           finite-difference schemes (-1: backward, 1: forward)*/
                                        void    *s,   /*!< Solver object of type #SolverContext */
                                        void    *m    /*!< MPI object of type #MPIContext */
                                     )
{
  SolverContext *solver = (SolverContext*) s;
  int           i, j, v;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  if ((!Df) || (!f)) {
    fprintf(stderr, "Error in FirstDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

#pragma omp parallel for schedule(auto) default(shared) private(i,j,v,index_outer,indexC)
  for (j=0; j<N_outer; j++) {
    _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    /* left boundary */
    for (i = -ghosts; i < -ghosts+1; i++) {
      int qC, qR, qRR;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR );
      indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qRR);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = 0.5 * (-3*f[qC*nvars+v]+4*f[qR*nvars+v]-f[qRR*nvars+v]);
    }
    /* interior */
    for (i = -ghosts+1; i < dim[dir]+ghosts-1; i++) {
      int qC, qL, qR;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = 0.5 * (f[qR*nvars+v]-f[qL*nvars+v]);
    }
    /* right boundary */
    for (i = dim[dir]+ghosts-1; i < dim[dir]+ghosts; i++) {
      int qLL, qL, qC;
      indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qLL);
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL );
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = 0.5 * (3*f[qC*nvars+v]-4*f[qL*nvars+v]+f[qLL*nvars+v]);
    }
  }
  
  return(0);
}
