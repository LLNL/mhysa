/*! @file ShallowWater2DEigen.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute left and right eigenvectors for the 2D shallow water equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <hypar.h>

/*! Compute the left eigenvections for the 2D shallow water equations. This function 
    just calls the macro #_ShallowWater2DLeftEigenvectors_ and is not used by any 
    functions within the 2D shallow water module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int ShallowWater2DLeftEigenvectors(
                            double *u, /*!< Conserved solution at a grid point */
                            double *L, /*!< Array of size nvar^2 = 3^2 to save the matrix of
                                            left eigenvectors in (row-major format). */
                            void   *p, /*!< Object of type #ShallowWater2D with physics-related variables */
                            int    dir /*!< Spatial dimension */
                           )
{
  ShallowWater2D *param  = (ShallowWater2D*)  p;
  _ShallowWater2DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

/*! Compute the right eigenvections for the 2D shallow water equations. This function 
    just calls the macro #_ShallowWater2DRightEigenvectors_ and is not used by any 
    functions within the 2D shallow water module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int ShallowWater2DRightEigenvectors(
                              double  *u, /*!< Conserved solution at a grid point */
                              double  *R, /*!< Array of size nvar^2 = 3^2 to save the matrix of 
                                               right eigenvectors in (row-major format). */
                              void    *p, /*!< Object of type #ShallowWater2D with physics-related variables */
                              int     dir /*!< Spatial dimension */
                            )
{
  ShallowWater2D *param  = (ShallowWater2D*)  p;
  _ShallowWater2DRightEigenvectors_(u,R,param,dir);
  return(0);
}
