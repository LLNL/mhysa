#include <stdlib.h>
#include <advectiondiffusionreaction.h>
#include <mpivars.h>
#include <hypar.h>

int LinearADRCleanup(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
 
  LinearADR *physics = solver->physics;
  if (physics->a) free(physics->a);

  return(0);
}
