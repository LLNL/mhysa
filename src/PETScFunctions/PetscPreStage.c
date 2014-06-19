#ifdef with_petsc

#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPreStage"

PetscErrorCode PetscPreStage(TS ts,PetscReal waqt)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#endif
