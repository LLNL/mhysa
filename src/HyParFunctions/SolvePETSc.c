/*! @file SolvePETSc.c
    @brief Integrate in time using PETSc
    Integrate the spatially discretized system in time using PETSc's TS module.\n
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html)
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "SolvePETSc"

/*! \brief Integrate in time with PETSc

    This function integrates the semi-discrete ODE (obtained from discretizing
    the PDE in space) using the time-integration module of PETSc 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html).
    The time-integration context is set up using the parameters specified in 
    the input file. However, they can also be specified using PETSc's command
    line inputs.\n
    \n
    See PETSc's documentation and examples for more details on how to use its
    TS module.
*/

int SolvePETSc(void *s, /*!< Solver object of type #HyPar */
               void *m  /*!< MPI object of type #MPIVariables */)
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  PetscErrorCode  ierr    = 0, d;
  TS              ts;     /* time integration object               */
  Vec             Y;      /* PETSc solution vector                 */
  Mat             A, B;   /* Jacobian and preconditioning matrices */
  TSType          time_scheme;  /* time integration method         */
  TSProblemType   ptype;  /* problem type - nonlinear or linear    */
  int             flag_mat_a = 0, flag_mat_b = 0;

  PetscFunctionBegin;

  /* Register custom time-integration methods, if specified */
  ierr = PetscRegisterTIMethods(mpi->rank); CHECKERR(ierr);
  if(!mpi->rank) printf("Setting up PETSc time integration... \n");

  /* create and set a PETSc context */
  PETScContext context;
  context.solver = solver;
  context.mpi    = mpi;
  context.tic    = 0;
  context.flag_hyperbolic     = _EXPLICIT_; 
  context.flag_hyperbolic_f   = _EXPLICIT_; 
  context.flag_hyperbolic_df  = _EXPLICIT_; 
  context.flag_parabolic      = _EXPLICIT_; 
  context.flag_source         = _EXPLICIT_; 
  context.flag_is_linear      = 0;

  /* create and initialize PETSc solution vector and other parameters */
  /* PETSc solution vector does not have ghost points */
  int total_size = 1;
  for (d=0; d<solver->ndims; d++) total_size *= (solver->dim_local[d]);
  total_size *= solver->nvars;
  ierr = VecCreate(MPI_COMM_WORLD,&Y); CHKERRQ(ierr);
  ierr = VecSetSizes(Y,total_size,PETSC_DECIDE); CHKERRQ(ierr);
  ierr = VecSetUp(Y); CHKERRQ(ierr);

  /* copy initial solution to PETSc's vector */
  ierr = TransferVecToPETSc(solver->u,Y,&context); CHECKERR(ierr);

  /* Define and initialize the time-integration object */
  ierr = TSCreate(MPI_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetDuration(ts,solver->n_iter,solver->dt*solver->n_iter); CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,solver->dt); CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  /* Define the right and left -hand side functions for each time-integration scheme */
  ierr = TSGetType(ts,&time_scheme); CHKERRQ(ierr);
  ierr = TSGetProblemType(ts,&ptype); CHKERRQ(ierr);
  if (!strcmp(time_scheme,TSARKIMEX)) {

    ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionIMEX,&context); CHKERRQ(ierr);
    ierr = TSSetIFunction  (ts,PETSC_NULL,PetscIFunctionIMEX,  &context); CHKERRQ(ierr);

    SNES     snes;
    KSP      ksp;
    PC       pc;
    SNESType snestype;
    ierr = TSGetSNES(ts,&snes); CHKERRQ(ierr);
    ierr = SNESGetType(snes,&snestype); CHKERRQ(ierr);

    /* Matrix-free representation of the Jacobian */
    flag_mat_a = 1;
    ierr = MatCreateShell(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,
                          PETSC_DETERMINE,&context,&A); CHKERRQ(ierr);
    if ((!strcmp(snestype,SNESKSPONLY)) || (ptype == TS_LINEAR)) {
      /* linear problem */
      context.flag_is_linear = 1;
      ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_Linear); CHKERRQ(ierr);
      ierr = SNESSetType(snes,SNESKSPONLY); CHKERRQ(ierr);
    } else {
      /* nonlinear problem */
      context.flag_is_linear = 0;
      context.jfnk_eps = 1e-7;
      ierr = PetscOptionsGetReal(NULL,"-jfnk_epsilon",&context.jfnk_eps,NULL); CHKERRQ(ierr);
      ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_JFNK); CHKERRQ(ierr);
    }
    ierr = MatSetUp(A); CHKERRQ(ierr);

    context.flag_use_precon = 0;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-with_pc",(PetscBool*)(&context.flag_use_precon),PETSC_NULL); CHKERRQ(ierr);

    if (context.flag_use_precon) {
      /* check if flux Jacobian of the physical model is defined */
      if (!solver->JFunction) {
        if (!mpi->rank) {
          fprintf(stderr,"Error in SolvePETSc(): solver->JFunction (point-wise flux Jacobian) must ");
          fprintf(stderr,"be defined for preconditioning.\n");
        }
        PetscFunctionReturn(1);
      }
      /* Set up preconditioner matrix */
      flag_mat_b = 1;
      ierr = MatCreateAIJ(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE,
                          (solver->ndims*2+1)*solver->nvars,NULL,
                          2*solver->ndims*solver->nvars,NULL,&B); CHKERRQ(ierr);
      ierr = MatSetBlockSize(B,solver->nvars);
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,B,PetscIJacobianIMEX,&context); CHKERRQ(ierr);
    } else {
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,A,PetscIJacobianIMEX,&context); CHKERRQ(ierr);
      /* Set PC (preconditioner) to none */
      ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
      ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
    }

    /* read the implicit/explicit flags for each of the terms for IMEX schemes */
    /* default -> hyperbolic - explicit, parabolic and source - implicit       */
    PetscBool flag = PETSC_FALSE;

    context.flag_hyperbolic     = _EXPLICIT_; 
    context.flag_hyperbolic_f   = _EXPLICIT_; 
    context.flag_hyperbolic_df  = _IMPLICIT_; 
    context.flag_parabolic      = _IMPLICIT_; 
    context.flag_source         = _IMPLICIT_; 

    if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {

      flag = PETSC_FALSE; 
      ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_f_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_f = _EXPLICIT_; 
      flag = PETSC_FALSE; 
      ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_f_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_f = _IMPLICIT_; 

      flag = PETSC_FALSE; 
      ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_df_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_df = _EXPLICIT_; 
      flag = PETSC_FALSE; 
      ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_df_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_df = _IMPLICIT_; 

    } else {

      flag = PETSC_FALSE; 
      ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
      if (flag == PETSC_TRUE) context.flag_hyperbolic = _EXPLICIT_; 
      flag = PETSC_FALSE; 
      ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
      if (flag == PETSC_TRUE) context.flag_hyperbolic = _IMPLICIT_; 

    }

    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-parabolic_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_parabolic = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-parabolic_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_parabolic = _IMPLICIT_; 

    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-source_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_source = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-source_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_source = _IMPLICIT_; 

    /* print out a summary of the treatment of each term */
    if (!mpi->rank) {
      printf("Implicit-Explicit time-integration:-\n");
      if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
        if (context.flag_hyperbolic_f == _EXPLICIT_)  printf("Hyperbolic (f-df) term: Explicit\n");
        else                                          printf("Hyperbolic (f-df) term: Implicit\n");
        if (context.flag_hyperbolic_df == _EXPLICIT_) printf("Hyperbolic (df)   term: Explicit\n");
        else                                          printf("Hyperbolic (df)   term: Implicit\n");
      } else {
        if (context.flag_hyperbolic == _EXPLICIT_)    printf("Hyperbolic        term: Explicit\n");
        else                                          printf("Hyperbolic        term: Implicit\n");
      }
      if (context.flag_parabolic == _EXPLICIT_)       printf("Parabolic         term: Explicit\n");
      else                                            printf("Parabolic         term: Implicit\n");
      if (context.flag_source    == _EXPLICIT_)       printf("Source            term: Explicit\n");
      else                                            printf("Source            term: Implicit\n");
    }

  } else if ((!strcmp(time_scheme,TSEULER)) || (!strcmp(time_scheme,TSRK)) || (!strcmp(time_scheme,TSSSP))) {
    
    ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionExpl,&context); CHKERRQ(ierr);

  } else {

    ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionImpl,&context); CHKERRQ(ierr);

    SNES     snes;
    KSP      ksp;
    PC       pc;
    SNESType snestype;
    ierr = TSGetSNES(ts,&snes); CHKERRQ(ierr);
    ierr = SNESGetType(snes,&snestype); CHKERRQ(ierr);

    /* Matrix-free representation of the Jacobian */
    flag_mat_a = 1;
    ierr = MatCreateShell(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,
                          PETSC_DETERMINE,&context,&A); CHKERRQ(ierr);
    if ((!strcmp(snestype,SNESKSPONLY)) || (ptype == TS_LINEAR)) {
      /* linear problem */
      context.flag_is_linear = 1;
      ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionImpl_Linear); CHKERRQ(ierr);
      ierr = SNESSetType(snes,SNESKSPONLY); CHKERRQ(ierr);
    } else {
      /* nonlinear problem */
      context.flag_is_linear = 0;
      context.jfnk_eps = 1e-7;
      ierr = PetscOptionsGetReal(NULL,"-jfnk_epsilon",&context.jfnk_eps,NULL); CHKERRQ(ierr);
      ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionImpl_JFNK); CHKERRQ(ierr);
    }
    ierr = MatSetUp(A); CHKERRQ(ierr);

    context.flag_use_precon = 0;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-with_pc",(PetscBool*)(&context.flag_use_precon),PETSC_NULL); CHKERRQ(ierr);

    if (context.flag_use_precon) {
      /* check if flux Jacobian of the physical model is defined */
      if (!solver->JFunction) {
        if (!mpi->rank) {
          fprintf(stderr,"Error in SolvePETSc(): solver->JFunction (point-wise flux Jacobian) must ");
          fprintf(stderr,"be defined for preconditioning.\n");
        }
        PetscFunctionReturn(1);
      }
      /* Set up preconditioner matrix */
      flag_mat_b = 1;
      ierr = MatCreateAIJ(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE,
                          (solver->ndims*2+1)*solver->nvars,NULL,
                          2*solver->ndims*solver->nvars,NULL,&B); CHKERRQ(ierr);
      ierr = MatSetBlockSize(B,solver->nvars);
      /* Set the RHSJacobian function for TS */
      ierr = TSSetRHSJacobian(ts,A,B,PetscRHSJacobian,&context); CHKERRQ(ierr);
    } else {
      /* Set the RHSJacobian function for TS */
      ierr = TSSetRHSJacobian(ts,A,A,PetscRHSJacobian,&context); CHKERRQ(ierr);
      /* Set PC (preconditioner) to none */
      ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
      ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
    }

  }

  /* Set pre/post-stage and post-timestep function */
  ierr = TSSetPreStep (ts,PetscPreTimeStep ); CHKERRQ(ierr);
  ierr = TSSetPreStage(ts,PetscPreStage    ); CHKERRQ(ierr);
  ierr = TSSetPostStage(ts,PetscPostStage  ); CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,PetscPostTimeStep); CHKERRQ(ierr);
  /* Set solution vector for TS */
  ierr = TSSetSolution(ts,Y); CHKERRQ(ierr);
  /* Set it all up */
  ierr = TSSetUp(ts); CHKERRQ(ierr);
  /* Set application context */
  ierr = TSSetApplicationContext(ts,&context); CHKERRQ(ierr);

  if (!mpi->rank) {
    if (context.flag_is_linear) printf("SolvePETSc(): Problem type is linear.\n");
    else                        printf("SolvePETSc(): Problem type is nonlinear.\n");
  }

  if (!mpi->rank) printf("** Starting PETSc time integration **\n");
  ierr = TSSolve(ts,Y); CHKERRQ(ierr);
  if (!mpi->rank) printf("** Completed PETSc time integration **\n");

  /* Get the number of time steps */
  ierr = TSGetTimeStepNumber(ts,&solver->n_iter); CHKERRQ(ierr);

  /* copy final solution from PETSc's vector */
  ierr = TransferVecFromPETSc(solver->u,Y,&context); CHECKERR(ierr);

  /* clean up */
  if (flag_mat_a) { ierr = MatDestroy(&A); CHKERRQ(ierr); }
  if (flag_mat_b) { ierr = MatDestroy(&B); CHKERRQ(ierr); }
  ierr = TSDestroy(&ts); CHKERRQ(ierr);
  ierr = VecDestroy(&Y); CHKERRQ(ierr);

  /* write a final solution file, if last iteration did not write one */
  if (context.tic) { 
    if (solver->PhysicsOutput) {
      IERR solver->PhysicsOutput(solver,mpi); CHECKERR(ierr);
    }
    IERR OutputSolution(solver,mpi); CHECKERR(ierr); 
  }
  /* calculate error if exact solution has been provided */
  IERR CalculateError(solver,mpi); CHECKERR(ierr);

  PetscFunctionReturn(0);
}

#endif
