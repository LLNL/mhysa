/*! @file Euler1DPreStep.c
    @author Debojyoti Ghosh
    @brief Contains the 1D Euler-specific function to be called at the beginning of each time step.
*/

#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! 
  1D Euler-specific function called at the beginning of each time-step
*/
int Euler1DPreStep(
                    double  *u,   /*!< Solution (conserved variables) */
                    void    *s,   /*!< Solver object of type #HyPar */
                    void    *m,   /*!< MPI object of type #MPIVariables */
                    double  waqt  /*!< Current solution time */
                  )
{
  return(0);
}
