/*! @file ShallowWater2DJacobian.c
    @author Debojyoti Ghosh
    @brief Contains the functions compute flux Jacobians for the 2D shallow water system
*/

#include <mathfunctions.h>
#include <physicalmodels/shallowwater2d.h>

/*! Function to compute the flux Jacobian of the 2D shallow water equations, given the
    solution at a grid point. The Jacobian is square matrix of size nvar=2, and 
    is returned as a 1D array (double) of 4 elements in row-major format.
*/
int ShallowWater2DJacobian(
                    double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 9 */
                    double  *u,   /*!< solution at a grid point (array of size nvar = 3) */
                    void    *p,   /*!< object containing the physics-related parameters */
                    int     dir,  /*!< spatial dimension (x/y) */
                    int     upw   /*!< 0 -> send back complete Jacobian, 
                                       1 -> send back Jacobian of right(+)-moving flux, 
                                      -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  ShallowWater2D  *param = (ShallowWater2D*) p;
  static double   R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                  L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _ShallowWater2DEigenvalues_      (u,D,param,dir);
  _ShallowWater2DLeftEigenvectors_ (u,L,param,dir);
  _ShallowWater2DRightEigenvectors_(u,R,param,dir);

  int aupw = absolute(upw), k;
  k = 0; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 4; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 8; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );

  MatMult3(_MODEL_NVARS_,DL,D,L);
  MatMult3(_MODEL_NVARS_,Jac,R,DL);

  return(0);
}
