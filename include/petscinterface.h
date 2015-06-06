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

  /* flags for Jacobian and preconditioning */
  int flag_jfnk_nopre;    /*!< use unpreconditioned Jacobian-free Newton-Krylov */
  int flag_jfnk_pre;      /*!< use preconditioned Jacobian-free Newton-Krylov   */

  /*! Flag to indicate if the system being solved for implicit time-integration is linear/nonlinear.
      It's decided on the basis of whether user specifies SNES type as ksponly.
  */
  int flag_is_linear;

} PETScContext;

/* Copy Functions */
/*! Function to copy the solution \a u in #HyPar to a Vec variable in PETSc 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/Vec.html) */
int TransferVecToPETSc    (double*,Vec,void*);
/*! Function to copy the solution \a u in #HyPar from a Vec variable in PETSc 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/Vec.html) */
int TransferVecFromPETSc  (double*,Vec,void*);
/*! Function to copy a matrix to a Mat variable in PETSc 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/Mat.html) */
int TransferMatToPETSc    (void*,Mat,void*);

/*! Function to register custom time-integration RK/ARKIMEX method, given the Butcher tableau. 
    \a Examples/PETScInputs has examples of input files needed to register a custom method. See
    \a Examples/README for more details. */
int PetscRegisterTIMethods (int);

/* Right and left -hand side functions */
/*! Compute the right-hand-side \a f(y) for explicit time integration of the ODE \a dy/dt = \a f(y).\n 
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSFunction.html 
*/
PetscErrorCode PetscRHSFunctionExpl (TS,PetscReal,Vec,Vec,void*);
/*! Compute the right-hand-side \a f(y) for implicit time integration of the ODE \a dy/dt = \a f(y).\n 
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSFunction.html 
*/
PetscErrorCode PetscRHSFunctionImpl (TS,PetscReal,Vec,Vec,void*);
/*! Compute the right-hand-side \a f(y) for implicit-explicit time integration of the ODE \a dy/dt - \a g(y) = \a f(y).\n 
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetRHSFunction.html 
*/
PetscErrorCode PetscRHSFunctionIMEX (TS,PetscReal,Vec,Vec,void*);
/*! Compute the left-hand-side \a dy/dt - \a g(y) for implicit-explicit time integration of the ODE \a dy/dt - \a g(y) = \a f(y).\n 
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIFunction.html 
*/
PetscErrorCode PetscIFunctionIMEX   (TS,PetscReal,Vec,Vec,Vec,void*);

/* Jacobian functions for left-hand side (IMEX) */
/*! Compute the Jacobian of the left-hand-side \a dy/dt - \a g(y) for the implicit-explicit time integration
    of the ODE \a dy/dt - \a g(y) = \a f(y), using the Jacobian-free Newtown-Krylov method with no
    preconditioning. \n\n**Note** that for the Jacobian-free Newton-Krylov approach, the Jacobian matrix is
    defined as the PETSc type MatShell (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
    This function just saves the shift and the current time. The action of the Jacobian matrix is defined through 
    #PetscJacobianFunctionIMEX_JFNK for nonlinear problems and #PetscJacobianFunctionIMEX_Linear for linear problems.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscIJacobianIMEX_JFNK_NoPre        (TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
/*! Compute the Jacobian of the left-hand-side \a dy/dt - \a g(y) for the implicit-explicit time integration
    of the ODE \a dy/dt - \a g(y) = \a f(y), using the preconditioned Jacobian-free Newtown-Krylov method. \n\n
    **Note** that for the Jacobian-free Newton-Krylov approach, the Jacobian matrix is
    defined as the PETSc type MatShell (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
    This function computes the preconditioning matrix and saves the shift and the current time. The action of the 
    Jacobian matrix is defined through #PetscJacobianFunctionIMEX_JFNK for nonlinear problems and 
    #PetscJacobianFunctionIMEX_Linear for linear problems.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscIJacobianIMEX_JFNK_Pre          (TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
/*! Compute the Jacobian and preconditioning matrix of the left-hand-side \a dy/dt - \a g(y) for 
    the implicit-explicit time integration of the ODE \a dy/dt - \a g(y) = \a f(y).\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscIJacobianIMEX_Jac_Pre           (TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
/*! Compute the Jacobian of the left-hand-side \a dy/dt - \a g(y) for the implicit-explicit time integration
    of the ODE \a dy/dt - \a g(y) = \a f(y), using the preconditioned Jacobian-free Newtown-Krylov method. \n\n
    **Note** that for the Jacobian-free Newton-Krylov approach, the Jacobian matrix is
    defined as the PETSc type MatShell (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
    The provided Jacobian matrix is used as the preconditioner and computed in this function. The shift and the current 
    time are saved. The action of the Jacobian matrix is defined through #PetscJacobianFunctionIMEX_JFNK for nonlinear 
    problems and #PetscJacobianFunctionIMEX_Linear for linear problems.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscIJacobianIMEX_JFNK_JacIsPre     (TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
/*! Compute the Jacobian of the left-hand-side \a dy/dt - \a g(y) for 
    the implicit-explicit time integration of the ODE \a dy/dt - \a g(y) = \a f(y). This is 
    used for the case where a separate preconditioning matrix is not defined.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscIJacobianIMEX_Jac_NoPre         (TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
/*! Function that defines the action of the Jacobian matrix on the input vector for a nonlinear problem
    using the Jacobian-free Newtown-Krylov approach.
    The Jacobian of the left-hand-side \a dy/dt - \a g(y) is computed for the implicit time-integration
    of the ODE \a dy/dt - \a g(y) = \a f(y).
    \sa #PetscIJacobianIMEX_JFNK_NoPre, #PetscIJacobianIMEX_JFNK_Pre
*/
PetscErrorCode PetscJacobianFunctionIMEX_JFNK       (Mat,Vec,Vec);             
/*! Function that defines the action of the Jacobian matrix on the input vector for a linear problem.
    The Jacobian of the left-hand-side \a dy/dt - \a g(y) is computed for the implicit time-integration
    of the ODE \a dy/dt - \a g(y) = \a f(y).
    \sa #PetscIJacobianIMEX_JFNK_NoPre, #PetscIJacobianIMEX_JFNK_Pre
*/
PetscErrorCode PetscJacobianFunctionIMEX_Linear     (Mat,Vec,Vec);             

/* Jacobian functions for right-hand side (implicit) */
/*! Compute the Jacobian of the right-hand-side \a f(y) for the implicit time integration
    of the ODE \a dy/dt = \a f(y), using the Jacobian-free Newtown-Krylov method with no
    preconditioning. \n\n**Note** that for the Jacobian-free Newton-Krylov approach, the Jacobian matrix is
    defined as the PETSc type MatShell (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
    This function just saves the shift and the current time. The action of the Jacobian matrix is defined through 
    #PetscJacobianFunctionIMEX_JFNK for nonlinear problems and #PetscJacobianFunctionIMEX_Linear for linear problems.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscRHSJacobianImpl_JFNK_NoPre      (TS,PetscReal,Vec,Mat,Mat,void*);
/*! Compute the Jacobian of the right-hand-side \a f(y) for the implicit time integration
    of the ODE \a dy/dt = \a f(y), using the preconditioned Jacobian-free Newtown-Krylov method. \n\n
    **Note** that for the Jacobian-free Newton-Krylov approach, the Jacobian matrix is
    defined as the PETSc type MatShell (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
    This function computes the preconditioning matrix and saves the shift and the current time. The action of the 
    Jacobian matrix is defined through #PetscJacobianFunctionIMEX_JFNK for nonlinear problems and 
    #PetscJacobianFunctionIMEX_Linear for linear problems.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscRHSJacobianImpl_JFNK_Pre        (TS,PetscReal,Vec,Mat,Mat,void*);
/*! Compute the Jacobian and preconditioning matrix of the right-hand-side \a f(y) for 
    the implicit time integration of the ODE \a dy/dt = \a f(y).\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscRHSJacobianImpl_Jac_Pre         (TS,PetscReal,Vec,Mat,Mat,void*);
/*! Compute the Jacobian of the right-hand-side \a f(y) for the implicit time integration
    of the ODE \a dy/dt = \a f(y), using the preconditioned Jacobian-free Newtown-Krylov method. \n\n
    **Note** that for the Jacobian-free Newton-Krylov approach, the Jacobian matrix is
    defined as the PETSc type MatShell (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html).
    The provided Jacobian matrix is used as the preconditioner and computed in this function. The shift and the current 
    time are saved. The action of the Jacobian matrix is defined through #PetscJacobianFunctionIMEX_JFNK for nonlinear 
    problems and #PetscJacobianFunctionIMEX_Linear for linear problems.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscRHSJacobianImpl_JFNK_JacIsPre   (TS,PetscReal,Vec,Mat,Mat,void*);
/*! Compute the Jacobian of the right-hand-side \a f(y) for 
    the implicit time integration of the ODE \a dy/dt = \a f(y). This is 
    used for the case where a separate preconditioning matrix is not defined.\n\n
    See: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSSetIJacobian.html
*/
PetscErrorCode PetscRHSJacobianImpl_Jac_NoPre       (TS,PetscReal,Vec,Mat,Mat,void*);
/*! Function that defines the action of the Jacobian matrix on the input vector for a nonlinear problem
    using the Jacobian-free Newtown-Krylov approach. The Jacobian of the right-hand-side \a f(y) is 
    computed for the implicit time-integration of the ODE \a dy/dt = \a f(y).
    \sa #PetscIJacobianImpl_JFNK_NoPre, #PetscIJacobianImpl_JFNK_Pre
*/
PetscErrorCode PetscJacobianFunctionImpl_JFNK       (Mat,Vec,Vec);             
/*! Function that defines the action of the Jacobian matrix on the input vector for a linear problem.
    The Jacobian of the right-hand-side \a f(y) is computed for the implicit time-integration of the 
    ODE \a dy/dt = \a f(y).
    \sa #PetscIJacobianImpl_JFNK_NoPre, #PetscIJacobianImpl_JFNK_Pre
*/
PetscErrorCode PetscJacobianFunctionImpl_Linear     (Mat,Vec,Vec);             

/* Other functions */
/*! Function called at the beginning of each stage computation in multi-stage time-integration methods. */
PetscErrorCode PetscPreStage        (TS,PetscReal);
/*! Function called at the end of each stage computation in multi-stage time-integration methods. */
PetscErrorCode PetscPostStage       (TS,PetscReal,PetscInt,Vec*);
/*! Function called at the beginning of each time step. */
PetscErrorCode PetscPreTimeStep     (TS);
/*! Function called at the end of each time step. */
PetscErrorCode PetscPostTimeStep    (TS);

