/*! @file LinearADRJacobian.c
    @author Debojyoti Ghosh
    @brief Function to compute the hyperbolic flux Jacobian for the linear advection-diffusion-reaction system.
*/

#include <mathfunctions.h>
#include <physicalmodels/linearadr.h>

/*! Function to compute the flux Jacobian for the hyperbolic (advection) part of the 
    linear-advection-diffusion-reaction model.
*/
int LinearADRJacobian(
                      double *Jac, /*!< Jacobian matrix of size 1 (nvar = 1) */
                      double *u,   /*!< solution at a grid point */
                      void *p,     /*!< object containing physics-related parameters */
                      int dir,     /*!< dimension (x/y/z) */
                      int upw      /*!< 0 -> send back complete Jacobian, 
                                        1 -> send back Jacobian of right(+)-moving flux, 
                                       -1 -> send back Jacobian of left(-)-moving flux*/
                     )
{
  LinearADR *param = (LinearADR*) p;
  *Jac =    (1-absolute(upw))*absolute(param->a[dir])
         +  absolute(upw) * (1+upw) * max(0,param->a[dir]) * 0.5
         -  absolute(upw) * (1-upw) * min(0,param->a[dir]) * 0.5 ;

  return(0);
}
