#include <stdlib.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

int Numa3DCleanup(void *s)
{
  HyPar  *solver  = (HyPar*)         s;
  Numa3D *physics = (Numa3D*)        solver->physics;

  free(physics->rho0);
  free(physics->P0);
  free(physics->T0);
  return(0);
}
