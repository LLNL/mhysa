#include <stdlib.h>
#include <physicalmodels/fppowersystem3bus.h>

int FPPowerSystem3BusCleanup(void *s)
{
  FPPowerSystem3Bus *physics = (FPPowerSystem3Bus*) s;
  free(physics->G);
  free(physics->B);
  free(physics->Ainv);
  return(0);
}
