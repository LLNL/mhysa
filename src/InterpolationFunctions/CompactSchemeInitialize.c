/*! @file CompactSchemeInitialize.c
    @brief Initializes the compact schemes
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
  This function initializes the compact finite-difference methods: allocates the arrays
  to store the tridiagonal system.
*/
int CompactSchemeInitialize(
                              void *s,      /*!< Solver object of type #HyPar */
                              void *m,      /*!< MPI object of type #MPIVariables */
                              char *type    /*!< Type of interpolation */
                           )
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  CompactScheme *compact   = (CompactScheme*) solver->compact;

  int nvars = solver->nvars;
  int ndims = solver->ndims;

  int size = 1, d;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+1);
  size *= solver->nvars;
  if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) size *= solver->nvars;

  compact->A = (double*) calloc (size, sizeof(double));
  compact->B = (double*) calloc (size, sizeof(double));
  compact->C = (double*) calloc (size, sizeof(double));
  compact->R = (double*) calloc (size, sizeof(double));

  compact->sendbuf = (double*) calloc (size, sizeof(double));
  compact->recvbuf = (double*) calloc (size, sizeof(double));

  return(0);
}
