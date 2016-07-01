/*! @file ApplyIBConditions.c
    @brief Apply immersed boundary conditions
    @author Debojyoti Ghosh
*/

#include <mpivars.h>
#include <hypar.h>

/*!
  Call the physics-specific function that applies the immersed boundary conditions
  on the immersed boundary points.
*/
int ApplyIBConditions(void    *s, /*!< Object of type #HyPar containing solver-related variables */
                      void    *m, /*!< Object of type #MPIVariables containing MPI-related variables */
                      double  *x, /*!< The solution vector to apply immersed BCs on */
                      double  waqt /*!< Current simulation time */
                     ) 
{
  HyPar           *solver   = (HyPar*)          s;
  MPIVariables    *mpi      = (MPIVariables*)   m;

  /* Apply immersed boundary conditions, if applicable */
  if (solver->flag_ib) IERR solver->IBFunction(solver,mpi,x,waqt);

  return(0);
}
