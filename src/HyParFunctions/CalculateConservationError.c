/*! @file CalculateConservationError.c
    @author Debojyoti Ghosh
    @brief Compute the conservation error.
*/
#include <math.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

/*! Calculates the error (L2) in conservation by computing the difference between
    the initial volume integral of the solution, and the sum of the current 
    volume integral and the time integral of the boundary flux from the start
    of the simulation to the current simulation time.
*/
int CalculateConservationError(
                                void *s, /*!< Solver object of type #HyPar */
                                void *m  /*!< MPI object of type #MPIVariables */
                              )
{
  HyPar         *solver = (HyPar*) s;
  int           v,nvars = solver->nvars;
  double        error;

  for (v=0; v<nvars; v++) {
    error =  (solver->VolumeIntegral[v]+solver->TotalBoundaryIntegral[v]-solver->VolumeIntegralInitial[v]) 
           * (solver->VolumeIntegral[v]+solver->TotalBoundaryIntegral[v]-solver->VolumeIntegralInitial[v]);
    solver->ConservationError[v] = sqrt(error);
  }
  
  return(0);
}
