/*! @file NavierStokes2DJacobian.c
    @author Debojyoti Ghosh
    @brief Contains the functions compute flux Jacobians for the 2D Navier-Stokes system
*/

#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>

/*! Function to compute the flux Jacobian of the 2D Navier-Stokes equations, given the
    solution at a grid point. The Jacobian is square matrix of size nvar=4, and 
    is returned as a 1D array (double) of 16 elements in row-major format.
*/
int NavierStokes2DJacobian(
                    double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 16 */
                    double  *u,   /*!< solution at a grid point (array of size nvar = 4) */
                    void    *p,   /*!< object containing the physics-related parameters */
                    int     dir,  /*!< dimension (0 -> x, 1 -> y) */
                    int     upw   /*!< 0 -> send back complete Jacobian, 
                                       1 -> send back Jacobian of right(+)-moving flux, 
                                      -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  NavierStokes2D *param = (NavierStokes2D*) p;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _NavierStokes2DEigenvalues_      (u,D,param,dir);
  _NavierStokes2DLeftEigenvectors_ (u,L,param,dir);
  _NavierStokes2DRightEigenvectors_(u,R,param,dir);

  int aupw = absolute(upw), k;
  k = 0;  D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 5;  D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 10; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 15; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );

  MatMult4(_MODEL_NVARS_,DL,D,L);
  MatMult4(_MODEL_NVARS_,Jac,R,DL);

  return(0);
}

/*! Function to compute the Jacobian of the fast flux (representing the acoustic waves) 
    of the 2D Navier-Stokes equations, given the solution at a grid point (see #_NavierStokes2DSetStiffFlux_, 
    #_NavierStokes2DSetStiffJac_). The Jacobian is square matrix of size nvar=4, and is returned as 
    a 1D array (double) of 16 elements in row-major format.
*/
int NavierStokes2DStiffJacobian(
                          double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 16 */
                          double  *u,   /*!< solution at a grid point (array of size nvar = 4) */
                          void    *p,   /*!< object containing the physics-related parameters */
                          int     dir,  /*!< dimension (0 -> x, 1 -> y) */
                          int     upw   /*!< 0 -> send back complete Jacobian, 
                                             1 -> send back Jacobian of right(+)-moving flux, 
                                            -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  NavierStokes2D *param = (NavierStokes2D*) p;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _NavierStokes2DEigenvalues_      (u,D,param,dir);
  _NavierStokes2DLeftEigenvectors_ (u,L,param,dir);
  _NavierStokes2DRightEigenvectors_(u,R,param,dir);

  int aupw = absolute(upw), k;
  k = 0;  D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 5;  D[k] = ( dir == _YDIR_ ? 0.0 : absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) ) );
  k = 10; D[k] = ( dir == _XDIR_ ? 0.0 : absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) ) );
  k = 15; D[k] = 0.0;

  MatMult4(_MODEL_NVARS_,DL,D,L);
  MatMult4(_MODEL_NVARS_,Jac,R,DL);

  return(0);
}
