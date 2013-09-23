#ifdef with_petsc

#include <petscinterface.h>

PetscErrorCode PetscPostStage(TS ts,PetscReal stagetime,PetscInt stageindex,Vec *Y)
{
  return(0);
}

#endif
