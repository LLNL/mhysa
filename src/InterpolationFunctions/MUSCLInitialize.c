/*! @file MUSCLInitialize.c
    @brief Initialize the 3rd order MUSCL scheme
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Initialize the 3rd order MUSCL scheme. The current implementation
    just needs to initialize the value of the parameter \f$\epsilon\f$,
    which is hard-coded to \f$10^{-3}\f$.
*/
int MUSCLInitialize(
                      void *s,  /*!< Solver object of type #HyPar */
                      void *m   /*!< MPI object of type #MPIVariables */
                   )
{
  MUSCLParameters  *muscl   = (MUSCLParameters*) s;
  /* hard coding these parameters for now */
  /* modify to read from an input file later */
  muscl->eps = 1e-3;
  return(0);
}
