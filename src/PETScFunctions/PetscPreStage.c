/*! @file PetscPreStage.c
    @brief Pre-time-integration-stage function
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPreStage"

/*! Function called before a stage in multi-stage time-integration methods */
PetscErrorCode PetscPreStage(
                              TS        ts,   /*!< Time integration object */
                              PetscReal waqt  /*!< Current simulation time */
                            )
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#endif
