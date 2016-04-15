/*! @file NavierStokes3DEigen.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute left and right eigenvectors for the 3D Navier Stokes equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>


/*! Compute the left eigenvections for the 3D Navier Stokes equations. This function 
    just calls the macro #_NavierStokes3DLeftEigenvectors_ and is not used by any 
    functions within the 3D Navier Stokes module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes3DLeftEigenvectors(
                                    double *u, /*!< Conserved solution at a grid point */
                                    double *L, /*!< Array of size nvar^2 = 5^2 to save the matrix of
                                                    left eigenvectors in (row-major format). */
                                    void   *p, /*!< Object of type #NavierStokes3D with physics-related variables */
                                    int    dir /*!< Spatial dimension (x, y, or z) */
                                  )
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  _NavierStokes3DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

/*! Compute the right eigenvections for the 3D Navier Stokes equations. This function 
    just calls the macro #_NavierStokes3DRightEigenvectors_ and is not used by any 
    functions within the 3D Navier Stokes module. However, it's necessary to define it 
    and provide it to the the solver object (#HyPar) so that it can then send it 
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes3DRightEigenvectors(
                                      double  *u, /*!< Conserved solution at a grid point */
                                      double  *R, /*!< Array of size nvar^2 = 5^2 to save the matrix of 
                                                       right eigenvectors in (row-major format). */
                                      void    *p, /*!< Object of type #NavierStokes3D with physics-related variables */
                                      int     dir /*!< Spatial dimension (x, y, or z) */
                                   )
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  _NavierStokes3DRightEigenvectors_(u,R,param,dir);
  return(0);
}
