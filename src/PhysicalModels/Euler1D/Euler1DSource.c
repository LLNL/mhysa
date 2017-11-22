/*! @file Euler1DSource.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the source terms for the 1D Euler equations.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the source terms for the 1D Euler equations. */
int Euler1DSource(
                  double  *source, /*!< Computed source terms (array size & layout same as u) */
                  double  *u,      /*!< Solution (conserved variables) */
                  void    *s,      /*!< Solver object of type #HyPar */
                  void    *m,      /*!< MPI object of type #MPIVariables */
                  double  t        /*!< Current solution time */
                 )
{
  HyPar         *solver = (HyPar* ) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  Euler1D       *param  = (Euler1D*) solver->physics;

  return(0);
}
