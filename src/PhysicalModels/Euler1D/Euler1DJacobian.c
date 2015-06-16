/*! @file Euler1DJacobian.c
    @author Debojyoti Ghosh
    @brief Contains the functions compute flux Jacobians for the 1D Euler system
*/

#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>

/*! Function to compute the flux Jacobian of the 1D Euler equations, given the
    solution at a grid point. The Jacobian is square matrix of size nvar=3, and 
    is returned as a 1D array (double) of 9 elements in row-major format.
*/
int Euler1DJacobian(
                    double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 9 */
                    double  *u,   /*!< solution at a grid point (array of size nvar = 3) */
                    void    *p,   /*!< object containing the physics-related parameters */
                    int     dir,  /*!< dimension (x/y/z) (not used, since this is 1D system) */
                    int     upw   /*!< 0 -> send back complete Jacobian, 
                                       1 -> send back Jacobian of right(+)-moving flux, 
                                      -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  Euler1D       *param = (Euler1D*) p;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _Euler1DEigenvalues_      (u,D,param,0);
  _Euler1DLeftEigenvectors_ (u,L,param,0);
  _Euler1DRightEigenvectors_(u,R,param,0);

  int aupw = absolute(upw), k;
  k = 0; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 4; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 8; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );

  MatMult3(3,DL,D,L);
  MatMult3(3,Jac,R,DL);

  return(0);
}

/*! Function to compute the Jacobian of the fast flux (representing the acoustic waves) 
    of the 1D Euler equations, given the solution at a grid point (see #_Euler1DSetStiffFlux_, 
    #_Euler1DSetStiffJac_). The Jacobian is square matrix of size nvar=3, and is returned as 
    a 1D array (double) of 9 elements in row-major format.
*/
int Euler1DStiffJacobian(
                          double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 9 */
                          double  *u,   /*!< solution at a grid point (array of size nvar = 3) */
                          void    *p,   /*!< object containing the physics-related parameters */
                          int     dir,  /*!< dimension (x/y/z) (not used, since this is 1D system) */
                          int     upw   /*!< 0 -> send back complete Jacobian, 
                                             1 -> send back Jacobian of right(+)-moving flux, 
                                            -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  Euler1D       *param = (Euler1D*) p;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], 
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _Euler1DEigenvalues_      (u,D,param,0);
  _Euler1DLeftEigenvectors_ (u,L,param,0);
  _Euler1DRightEigenvectors_(u,R,param,0);

  int aupw = absolute(upw), k;
  k = 0; D[k] = 0.0; /* remove the entropy eigenmode */
  k = 4; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 8; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );

  MatMult3(3,DL,D,L);
  MatMult3(3,Jac,R,DL);

  return(0);
}
