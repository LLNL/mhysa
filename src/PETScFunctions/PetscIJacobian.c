/*! @file PetscIJacobian.c
    @brief Contains the functions required for Jacobian computations for implicit time-integration
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscIJacobian"
/*! 
    Compute the Jacobian for implicit time integration of the governing equations: 
    The ODE, obtained after discretizing the governing PDE in space, is expressed as follows:
    \f{equation}{
      \frac {d{\bf U}}{dt} = {\bf F}\left({\bf U}\right) 
      \Rightarrow \dot{\bf U} - {\bf F}\left({\bf U}\right) = 0, 
    \f}
    where \f${\bf F}\f$ is the spatially discretized right-hand-side, and \f${\bf U}\f$ 
    represents the entire solution vector. The Jacobian is thus given by:
    \f{equation}{
      {\bf J} = \left[\alpha{\bf I} - \frac {\partial {\bf F}} {\partial {\bf U}} \right]
    \f}
    where \f$\alpha\f$ is the shift coefficient (#PETScContext::shift) of the time integration method.

    \b Matrix-free \b representation: The Jacobian is computed using a matrix-free approach, where the 
    entire Jacobian is never assembled, computed, or stored. Instead, the action of the Jacobian on
    a vector is defined through functions (PetscJacobianFunction_JFNK() for nonlinear systems, and
    PetscJacobianFunction_Linear() for linear systems). This approach works well with iterative 
    solvers. Thus, this function just does the following:
    + Saves the #PETScContext::shift (\f$\alpha\f$) and #PETScContext::waqt (current simulation time)
      to the application context (so that PetscJacobianFunction_JFNK() or PetscJacobianFunction_Linear() 
      can access these values).
    + If a preconditioner is being used, calls the function to compute the preconditioning matrix.

    \b Notes:
    + The Jacobian is defined as the PETSc type MatShell 
      (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSHELL.html)
    + \a Y and \a Ydot in the code are \f${\bf U}\f$ and \f$\dot{\bf U}\f$, respectively. PETsc denotes
      the state vector with \f${\bf Y}\f$ in its time integrators.
    + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html for documentation on
      PETSc's time integrators.
    + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
      the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
      or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscIJacobian(
                               TS ts,        /*!< Time stepping object (see PETSc TS)*/
                               PetscReal t,  /*!< Current time */
                               Vec Y,        /*!< Solution vector */
                               Vec Ydot,     /*!< Time-derivative of solution vector */
                               PetscReal a,  /*!< Shift */
                               Mat A,        /*!< Jacobian matrix */
                               Mat B,        /*!< Preconditioning matrix */
                               void *ctxt    /*!< Application context */
                             )
{
  PETScContext *context = (PETScContext*) ctxt;
  HyPar        *solver  = context->solver;
  _DECLARE_IERR_;

  PetscFunctionBegin;
  solver->count_IJacobian++;
  context->shift = a;
  context->waqt  = t;
  /* Construct preconditioning matrix */
  if (context->flag_use_precon) { IERR PetscComputePreconMatImpl(B,Y,context); CHECKERR(ierr); }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunction_JFNK"
/*! 
    Computes the action of the Jacobian on a vector: See documentation for PetscIJacobian()
    for the definition of the Jacobian. This function computes its action on a vector
    for a nonlinear system by taking the directional derivative, as follows:
    \f{equation}{
      {\bf f} = {\bf J}\left({\bf U}_0\right){\bf y} = \frac {\partial \mathcal{F}\left({\bf U}\right)} {\partial {\bf U}}\left({\bf U}_0\right) {\bf y} \approx \frac{1}{\epsilon} \left[ \mathcal{F}\left({\bf U}_0+\epsilon{\bf y}\right)-\mathcal{F}\left({\bf U}_0\right) \right]
    \f}
    In the context of the implicit time integration of the ODE given by
    \f{equation}{
      \frac {d {\bf U}} {dt} = {\bf F}\left({\bf U}\right) \Rightarrow \frac {d {\bf U}} {dt} - {\bf F}\left({\bf U}\right) = 0,
    \f}
    we have that
    \f{equation}{
      \mathcal{F}\left({\bf U}\right) \equiv \dot{\bf U} - {\bf F}\left({\bf U}\right)
      \Rightarrow {\bf J}\left({\bf U}\right) = \left[\alpha{\bf I} - \frac {\partial {\bf F}\left({\bf U}\right)} {\partial {\bf U}} \right].
    \f}
    where \f$\alpha\f$ (#PETScContext::shift) is the shift variable specific to the time integration method.
    So this function computes
    \f{equation}{
      {\bf f} = \alpha {\bf y} - \frac{1}{\epsilon} \left[ {\bf F}\left({\bf U}_0+\epsilon{\bf y}\right)-{\bf F}\left({\bf U}_0\right) \right]
    \f}
    In the code, \f${\bf y}\f$ is \a Y, \f$\bf f\f$ is \a F, and \f${\bf U}_0\f$ is \a #HyPar::uref (the reference solution at which the 
    nonlinear Jacobian is computed). See papers on Jacobian-free Newton-Krylov (JFNK) methods to understand how \f$\epsilon\f$ is computed.

    \b Notes:
    + For a nonlinear spatial discretization, the right-hand-side \b must be computed without the nonlinearity
      (i.e. with a previously computed or "frozen" discretization operator). This ensures that the Jacobian being
      computed is consistent.
    + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html for documentation on
      PETSc's time integrators.
    + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
      the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
      or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscJacobianFunction_JFNK(
                                           Mat Jacobian, /*!< Jacobian matrix */
                                           Vec Y,        /*!< Input vector */
                                           Vec F         /*!< Output vector (Jacobian times input vector) */
                                         )
{
  PETScContext    *context = NULL;
  HyPar           *solver  = NULL;
  MPIVariables    *mpi     = NULL;
  int             ierr     = 0, d;

  PetscFunctionBegin;

  ierr   = MatShellGetContext(Jacobian,&context); CHKERRQ(ierr);
  solver = context->solver;
  mpi    = context->mpi;
  solver->count_IJacFunction++;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u       = solver->u;
  double *uref    = solver->uref;
  double *rhsref  = solver->rhsref;
  double *rhs     = solver->rhs;

  double t = context->waqt; /* current stage/step time */

  double normY;
  ierr = VecNorm(Y,NORM_2,&normY);                                        CHKERRQ(ierr);

  if (normY < 1e-16) {

    /* F = 0 */
    ierr = VecZeroEntries(F);                                             CHKERRQ(ierr);
    /* [J]Y = aY - F(Y) */
    ierr = VecAXPBY(F,context->shift,0,Y);                                CHKERRQ(ierr);

  } else {
    
    double epsilon =  context->jfnk_eps / normY;
    /* copy solution from PETSc vector */
    ierr = TransferVecFromPETSc(u,Y,context);                             CHECKERR(ierr);
    _ArrayAYPX_(uref,epsilon,u,size*solver->nvars);
    /* apply boundary conditions and exchange data over MPI interfaces */
    ierr = solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);          CHECKERR(ierr);
    ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                   solver->ghosts,mpi,u);                 CHECKERR(ierr);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);  CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferVecToPETSc(rhs,F,context);               CHECKERR(ierr);
    /* [J]Y = aY - F(Y) */
    ierr = VecAXPBY(F,context->shift,(-1.0/epsilon),Y);     CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunction_Linear"
/*! 
    Computes the action of the Jacobian on a vector: See documentation for PetscIJacobian()
    for the definition of the Jacobian. This function computes its action on a vector
    for a linear system by taking the directional derivative, as follows:
    \f{equation}{
      {\bf f} = {\bf J}{\bf y} = \frac {\partial \mathcal{F}\left({\bf U}\right)} {\partial {\bf U}} {\bf y} = \left[ \mathcal{F}\left({\bf U}_0+{\bf y}\right)-\mathcal{F}\left({\bf U}_0\right) \right],
    \f}
    where \f$\mathcal{F}\f$ is linear, and thus \f${\bf J}\f$ is a constant (\f$\mathcal{F}\left({\bf y}\right) = {\bf J}{\bf y}\f$). 
    In the context of the implicit time integration of a linear 
    ODE given by
    \f{equation}{
      \frac {d {\bf U}} {dt} = {\bf F}\left({\bf U}\right) \Rightarrow \frac {d {\bf U}} {dt} - {\bf F}\left({\bf U}\right) = 0,
    \f}
    we have that
    \f{equation}{
      \mathcal{F}\left({\bf U}\right) \equiv \dot{\bf U} - {\bf F}\left({\bf U}\right)
      \Rightarrow {\bf J}\left({\bf U}\right) = \left[\alpha{\bf I} - \frac {\partial {\bf F}\left({\bf U}\right)} {\partial {\bf U}} \right].
    \f}
    where \f$\alpha\f$ (#PETScContext::shift) is the shift variable specific to the time integration method.
    So this function computes
    \f{equation}{
      {\bf f} = \alpha {\bf y} - \left[ {\bf F}\left({\bf U}_0+{\bf y}\right)-{\bf F}\left({\bf U}_0\right) \right]
    \f}
    In the code, \f${\bf y}\f$ is \a Y, \f$\bf f\f$ is \a F, and \f${\bf U}_0\f$ is \a #HyPar::uref (the reference solution at which the 
    nonlinear Jacobian is computed).

    Since \f$\mathcal{F}\f$ is linear, 
    \f{equation}{
      {\bf J}{\bf y} = \left[ \mathcal{F}\left({\bf U}_0+{\bf y}\right)-\mathcal{F}\left({\bf U}_0\right) \right] 
                     = \mathcal{F}\left({\bf y}\right).
    \f}
    However, the Jacobian is not computed as \f$\mathcal{F}\left({\bf y}\right)\f$ because of the following reason:
    this function is used by the equation solver within the implicit time integrator in PETSc, and \f${\bf y}\f$ 
    represents the change in the solution, i.e. \f$\Delta {\bf U}\f$, and not the solution \f$\bf U\f$. Thus, 
    evaluating \f$\mathcal{F}\left({\bf y}\right)\f$ using #HyPar::HyperbolicFunction, #HyPar::ParabolicFunction,
    and #HyPar::SourceFunction is incorrect since these functions expect \f$\bf U\f$ as the input.

    \b Notes:
    + For a nonlinear spatial discretization, the right-hand-side \b must be computed without the nonlinearity
      (i.e. with a previously computed or "frozen" discretization operator). This ensures that the Jacobian being
      computed is consistent, and is truly linear.
    + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html for documentation on
      PETSc's time integrators.
    + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
      the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
      or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscJacobianFunction_Linear(
                                             Mat Jacobian, /*!< Jacobian matrix */
                                             Vec Y,        /*!< Input vector */
                                             Vec F         /*!< Output vector (Jacobian times input vector */
                                           )
{
  PETScContext    *context = NULL;
  HyPar           *solver  = NULL;
  MPIVariables    *mpi     = NULL;
  int             ierr     = 0, d;

  PetscFunctionBegin;

  ierr   = MatShellGetContext(Jacobian,&context); CHKERRQ(ierr);
  solver = context->solver;
  mpi    = context->mpi;
  solver->count_IJacFunction++;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  double *u       = solver->u;
  double *uref    = solver->uref;
  double *rhsref  = solver->rhsref;
  double *rhs     = solver->rhs;

  double t = context->waqt; /* current stage/step time */

  double normY;
  ierr = VecNorm(Y,NORM_2,&normY);                                        CHKERRQ(ierr);

  if (normY < 1e-16) {

    /* F = 0 */
    ierr = VecZeroEntries(F);                                             CHKERRQ(ierr);
    /* [J]Y = aY - F(Y) */
    ierr = VecAXPBY(F,context->shift,0,Y);                                CHKERRQ(ierr);

  } else {
    
    /* copy solution from PETSc vector */
    ierr = TransferVecFromPETSc(u,Y,context);                             CHECKERR(ierr);
    _ArrayAYPX_(uref,1.0,u,size*solver->nvars);
    /* apply boundary conditions and exchange data over MPI interfaces */
    ierr = solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);          CHECKERR(ierr);
    ierr = MPIExchangeBoundariesnD(solver->ndims,solver->nvars,solver->dim_local,
                                   solver->ghosts,mpi,u);                 CHECKERR(ierr);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    ierr = solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);  CHECKERR(ierr);
    _ArrayAXPY_(solver->hyp,-1.0,rhs,size*solver->nvars);
    ierr = solver->ParabolicFunction (solver->par,u,solver,mpi,t);        CHECKERR(ierr);
    _ArrayAXPY_(solver->par, 1.0,rhs,size*solver->nvars);
    ierr = solver->SourceFunction    (solver->source,u,solver,mpi,t);     CHECKERR(ierr);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->nvars);
    /* Transfer RHS to PETSc vector */
    ierr = TransferVecToPETSc(rhs,F,context);     CHECKERR(ierr);
    /* [J]Y = aY - F(Y) */
    ierr = VecAXPBY(F,context->shift,-1.0,Y);     CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#endif
