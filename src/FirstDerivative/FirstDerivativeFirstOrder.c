/*! @file FirstDerivativeFirstOrder.c
    @author Debojyoti Ghosh
    @brief First order approximation to the first derivative
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <firstderivative.h>

#include <mpivars.h>
#include <hypar.h>
typedef MPIVariables  MPIContext;
typedef HyPar         SolverContext;

#ifdef with_omp
#include <omp.h>
#endif

/*! Computes the first-order finite-difference approximation to the first derivative 
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} f_{i+1} - f_i  & {\rm bias} = 1 \\ f_i - f_{i-1}  & {\rm bias} = -1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + The first derivative is computed at the ghost points too. Thus, biased schemes are used
      at and near the boundaries.
    + \b Df and \b f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.
*/
int FirstDerivativeFirstOrder(
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
      int qC, qR;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR );
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = f[qR*nvars+v]-f[qC*nvars+v];
    }
    /* interior */
    for (i = -ghosts+1; i < dim[dir]+ghosts-1; i++) {
      int qC, qL, qR;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = max(bias,0)*f[qR*nvars+v]-bias*f[qC*nvars+v]+min(bias,0)*f[qL*nvars+v];
    }
    /* right boundary */
    for (i = dim[dir]+ghosts-1; i < dim[dir]+ghosts; i++) {
      int qL, qC;
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL );
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = f[qC*nvars+v]-f[qL*nvars+v];
    }
  }
  
  return(0);
}
