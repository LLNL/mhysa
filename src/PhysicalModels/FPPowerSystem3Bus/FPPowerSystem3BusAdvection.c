/*! @file FPPowerSystem3BusAdvection.c
    @author Debojyoti Ghosh
    @brief Function to compute the advection term for the #FPPowerSystem3Bus system
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

/*! Compute the advection term for the #FPPowerSystem3Bus system: Since the advection
    coefficient is a function of x and not the solution, here the flux is set to the
    solution. The advection velocity is multiplied during upwinding FPPowerSystem3BusUpwind().
*/
int FPPowerSystem3BusAdvection(
                                double  *f,   /*!< Array to hold the computed flux vector (same layout as u) */
                                double  *u,   /*!< Array with the solution vector */
                                int     dir,  /*!< Spatial dimension for which to compute the flux */
                                void    *s,   /*!< Solver object of type #HyPar */
                                double  t     /*!< Current simulation time */
                              )
{
  HyPar *solver = (HyPar*) s;
  _ArrayCopy1D_(u,f,solver->npoints_local_wghosts);
  return(0);
}
