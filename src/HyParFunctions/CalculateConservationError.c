#include <math.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

int CalculateConservationError(void *s,void *m)
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           v,nvars = solver->nvars;
  double        error;

  for (v=0; v<nvars; v++) {
    error =  (solver->VolumeIntegral[v]+solver->TotalBoundaryIntegral[v]-solver->VolumeIntegralInitial[v]) 
           * (solver->VolumeIntegral[v]+solver->TotalBoundaryIntegral[v]-solver->VolumeIntegralInitial[v]);
    solver->ConservationError[v] = sqrt(error);
  }
  
  return(0);
}
