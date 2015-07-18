/*! @file ExactSolution.c
    @author Debojyoti Ghosh
    @brief Read in exact solution, if available.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <io.h>
#include <mpivars.h>
#include <hypar.h>

/*! Read in the exact solution, if available. */
int ExactSolution(
                    void    *s,     /*!< Solver object of type #HyPar */
                    void    *m,     /*!< MPI object of type #MPIVariables */
                    double  *uex,   /*!< Array to hold the exact solution, if available */
                    int     *flag   /*!< Flag to indicate if exact solution was available */
                 )
{
  HyPar  *solver = (HyPar*) s;
  ReadArray(solver->ndims,solver->nvars,solver->dim_global,solver->dim_local,
            solver->ghosts,solver,m,NULL,uex,"exact",flag);
  return(0);
}
