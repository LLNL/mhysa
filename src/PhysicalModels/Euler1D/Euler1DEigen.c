#if 0

/*! @file Euler1DEigen.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute left and right eigenvectors for the 1D Euler equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! Compute the left eigenvections for the 1D Euler equations. This function 
    just calls the macro #_Euler1DLeftEigenvectors_ and is not used by any 
    functions within the 1D Euler module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int Euler1DLeftEigenvectors(
                            double *u, /*!< Conserved solution at a grid point */
                            double *L, /*!< Array of size nvar^2 to save the matrix of
                                            left eigenvectors in (row-major format). */
                            void   *p, /*!< Object of type #Euler1D with physics-related variables */
                            int    dir /*!< Spatial dimension (not used, since this is a 1D system) */
                           )
{
  Euler1D *param  = (Euler1D*)  p;
  _Euler1DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

/*! Compute the right eigenvections for the 1D Euler equations. This function 
    just calls the macro #_Euler1DRightEigenvectors_ and is not used by any 
    functions within the 1D Euler module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int Euler1DRightEigenvectors(
                              double  *u, /*!< Conserved solution at a grid point */
                              double  *R, /*!< Array of size nvar^2 to save the matrix of 
                                               right eigenvectors in (row-major format). */
                              void    *p, /*!< Object of type #Euler1D with physics-related variables */
                              int     dir /*!< Spatial dimension (not used, since this is a 1D system) */
                            )
{
  Euler1D *param  = (Euler1D*)  p;
  _Euler1DRightEigenvectors_(u,R,param,dir);
  return(0);
}

#endif
