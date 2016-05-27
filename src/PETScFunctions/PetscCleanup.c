/*! @file PetscCleanup.c
    @author Debojyoti Ghosh
    @brief Cleans up allocations in PETSc interface.
*/

#ifdef with_petsc

#include <stdlib.h>
#include <petscinterface.h>

/*! Clean up allocations in the PETSc interface */
int PetscCleanup(void *obj /*!< Object of type #PETScContext */)
{
  PETScContext *ctxt = (PETScContext*) obj;
  if (ctxt->globalDOF) free(ctxt->globalDOF);
  return(0);
}

#endif
