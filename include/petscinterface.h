/* include PETSc header files */
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscts.h>

#define _MAX_STRING_SIZE_ 500

typedef struct _petsccontext_ {
  void *solver;
  void *mpi;

  PetscReal shift;
} PETScContext;

/* Copy Functions */
int TransferToPETSc    (double*,Vec,void*);
int TransferFromPETSc  (double*,Vec,void*);

/* Register custom time-integration RK/ARKIMEX method */
// int PetscRegisterTIMethods (char*,int);

/* Right and left -hand side functions */
PetscErrorCode PetscRHSFunctionExpl (TS,PetscReal,Vec,Vec,void*);
//PetscErrorCode PetscRHSFunctionIMEX (TS,PetscReal,Vec,Vec,void*);
//PetscErrorCode PetscIFunctionIMEX   (TS,PetscReal,Vec,Vec,Vec,void*);

/* Jacobian functions for left-hand side */
//PetscErrorCode PetscIJacobianIMEX         (TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);                            
//PetscErrorCode PetscJacobianFunctionIMEX  (Mat,Vec,Vec);             

/* Other functions */
PetscErrorCode PetscPreStage        (TS,PetscReal);
PetscErrorCode PetscPostStage       (TS,PetscReal,PetscInt,Vec*);
PetscErrorCode PetscPostTimeStep    (TS);

