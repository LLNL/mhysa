/*! @file petscinterface.h
    @brief Contains variables and function definitions for time integration with PETSc (http://www.mcs.anl.gov/petsc/)
    @author Debojyoti Ghosh
 */

/* include PETSc header files */
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscts.h>

/* some definitions */

/*! Maximum size of character strings */
#define _MAX_STRING_SIZE_ 500
/*! Explicit time integration */
#define _EXPLICIT_  0
/*! Implicit time integration */
#define _IMPLICIT_  1

/*! \def PETScContext
    \brief Structure containing the variables for time-integration with PETSc
 * This structure contains all the variables needed to integration in time using
 * PETSc's TS module (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html).
*/

/*! \brief  Structure containing the variables for time-integration with PETSc
 *
 * This structure contains all the variables needed to integration in time using
 * PETSc's TS module (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html).
*/
typedef struct _petsccontext_ {
  /*! object of type #HyPar containing the solver context */
  void *solver;
  /*! object of type #MPIVariables containing the MPI context */
  void *mpi;

  /*! The shift variable in implicit time-integration */
  PetscReal shift;
  /*! Current time */
  double    waqt;
  /*! A counter variable */
  int       tic;

  /*! Number of computational points, i.e., total number of grid points
      not counting the ghost points and blanked out points (for eg, inside
      immersed bodies) */
  int npoints;
  /*! index list of points whose size is (#HyPar::ndims+1)*#PETScContext::npoints.
      For each point, it stores its ndims-dimensional index and its 1D index 
      of its location in the array #HyPar::u.
  */
  int *points;

  /* flags for implicit treatment */
  /*! Flag to indicate if hyperbolic term is treated implicitly or explicitly    */
  int flag_hyperbolic;    
  /*! Flag to indicate if the split hyperbolic term (\a f - \a df) is treated implicitly or explicitly  */
  int flag_hyperbolic_f;  
  /*! Flag to indicate if the split hyperbolic term \a df is treated implicitly or explicitly  */
  int flag_hyperbolic_df; 
  /*! Flag to indicate if the parabolic term is treated implicitly or explicitly      */
  int flag_parabolic;     
  /*! Flag to indicate if the source term is treated implicitly or explicitly         */
  int flag_source;        

  /*! Flag to turn on preconditioning (off by default) */
  int flag_use_precon;

  /*! Flag to indicate if the system being solved for implicit time-integration is linear/nonlinear. */
  int flag_is_linear;

  /*! \f$\epsilon\f$ parameter for the Jacobian-free Newton-Krylov directional derivative computation */
  double jfnk_eps;

  /*! An essentialy integer array of the same size and layout at the solution 
      (with ghost points) containing the global DOF index for each grid point. 
      It is declared as a \a double type so that its calculation can use functions
      defined for double-type variables.
      \sa PetscGlobalDOF() */
  double *globalDOF;

} PETScContext;

/* Copy Functions */
int TransferVecToPETSc    (double*,Vec,void*);
int TransferVecFromPETSc  (double*,Vec,void*);
int TransferMatToPETSc    (void*,Mat,void*);

int PetscRegisterTIMethods (int);

/* Right and left -hand side functions */
PetscErrorCode PetscRHSFunctionExpl (TS,PetscReal,Vec,Vec,void*);
PetscErrorCode PetscRHSFunctionImpl (TS,PetscReal,Vec,Vec,void*);
PetscErrorCode PetscRHSFunctionIMEX (TS,PetscReal,Vec,Vec,void*);
PetscErrorCode PetscIFunctionIMEX   (TS,PetscReal,Vec,Vec,Vec,void*);

PetscErrorCode PetscIJacobianIMEX(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
PetscErrorCode PetscJacobianFunctionIMEX_JFNK       (Mat,Vec,Vec);             
PetscErrorCode PetscJacobianFunctionIMEX_Linear     (Mat,Vec,Vec);

PetscErrorCode PetscIJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
PetscErrorCode PetscJacobianFunction_JFNK  (Mat,Vec,Vec);             
PetscErrorCode PetscJacobianFunction_Linear(Mat,Vec,Vec);

int PetscComputePreconMatIMEX(Mat,Vec,void*);
int PetscComputePreconMatImpl(Mat,Vec,void*);

int PetscGlobalDOF(void*);
int PetscCleanup(void*);
int PetscCreatePointList(void*);

/* Other functions */
PetscErrorCode PetscPreStage        (TS,PetscReal);
PetscErrorCode PetscPostStage       (TS,PetscReal,PetscInt,Vec*);
PetscErrorCode PetscPreTimeStep     (TS);
PetscErrorCode PetscPostTimeStep    (TS);

/*! Function to compute any error estimates, if available */
PetscErrorCode PetscTimeError       (TS);

