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

int SolvePETSc(void *s,void *m)
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  PetscErrorCode  ierr    = 0, d;
  TS              ts;     /* time integration object               */
  Vec             Y;      /* PETSc solution vector                 */
  Mat             A, B;   /* Jacobian and preconditioning matrices */
  TSType          time_scheme;
  int             flag_mat_b = 0;

  PetscFunctionBegin;

  /* Register custom time-integration methods, if specified */
  ierr = PetscRegisterTIMethods(mpi->rank);                               CHECKERR(ierr);
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

  /* create and initialize PETSc solution vector and other parameters */
  /* PETSc solution vector does not have ghost points */
  int total_size = 1;
  for (d=0; d<solver->ndims; d++) total_size *= (solver->dim_local[d]);
  total_size *= solver->nvars;
  ierr = VecCreate(MPI_COMM_WORLD,&Y);                                    CHKERRQ(ierr);
  ierr = VecSetSizes(Y,total_size,PETSC_DECIDE);                          CHKERRQ(ierr);
  ierr = VecSetUp(Y);                                                     CHKERRQ(ierr);

  /* copy initial solution to PETSc's vector */
  ierr = TransferVecToPETSc(solver->u,Y,&context);                        CHECKERR(ierr);

  /* Define and initialize the time-integration object */
  ierr = TSCreate(MPI_COMM_WORLD,&ts);                                    CHKERRQ(ierr);
  ierr = TSSetDuration(ts,solver->n_iter,solver->dt*solver->n_iter);      CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,solver->dt);                         CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);             CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);                                            CHKERRQ(ierr);

  /* Define the right and left -hand side functions for each time-integration scheme */
  ierr = TSGetType(ts,&time_scheme);                                      CHKERRQ(ierr);
  if (   (!strcmp(time_scheme,TSEULER))
      || (!strcmp(time_scheme,TSSSP  ))
      || (!strcmp(time_scheme,TSRK   )) ){
    
    ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionExpl,&context); CHKERRQ(ierr);

  } else if (!strcmp(time_scheme,TSARKIMEX)) {

    ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionIMEX,&context); CHKERRQ(ierr);
    ierr = TSSetIFunction  (ts,PETSC_NULL,PetscIFunctionIMEX,  &context); CHKERRQ(ierr);

    /* read in the Jacobian-evaluation related flags */
    context.flag_jfnk_nopre = 0;
    context.flag_jfnk_pre   = 0;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-jfnk_nopre",(PetscBool*)(&context.flag_jfnk_nopre),PETSC_NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(PETSC_NULL,"-jfnk_pre"  ,(PetscBool*)(&context.flag_jfnk_pre)  ,PETSC_NULL); CHKERRQ(ierr);

    if ( ((!solver->JFunction) && (!solver->PFunction)) || (context.flag_jfnk_nopre) ){

      /* Physical model does not specify Jacobian and preconditioning functions, or
         user input flag specifies:
         Use: Unpreconditioned Jacobian-free Newton-Krylov
      */

      if (!mpi->rank) {
        printf("No Jacobian or preconditioner provided. ");
        printf("Using the unpreconditioned Jacobian-free Newton-Krylov approach.\n");
      }
      /* set pre-conditioner to none for MatShell */
      SNES     snes;
      KSP      ksp;
      PC       pc;
      SNESType snestype;
      ierr = TSGetSNES(ts,&snes);                                                   CHKERRQ(ierr);

      ierr = SNESGetType(snes,&snestype);                                           CHKERRQ(ierr);
      ierr = MatCreateShell(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,
                            PETSC_DETERMINE,&context,&A);                           CHKERRQ(ierr);
      if (!strcmp(snestype,SNESKSPONLY)) {
        context.flag_is_linear = 1;
        ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_Linear);
                                                                                    CHKERRQ(ierr);
      } else {
        context.flag_is_linear = 0;
        ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_JFNK);
                                                                                    CHKERRQ(ierr);
      }
      ierr = MatSetUp(A);                                                           CHKERRQ(ierr);
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,A,PetscIJacobianIMEX_JFNK_NoPre,&context);         CHKERRQ(ierr);
      /* Set PC (preconditioner) to none */
      ierr = SNESGetKSP(snes,&ksp);                                                 CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);                                                     CHKERRQ(ierr);
      ierr = PCSetType(pc,PCNONE);                                                  CHKERRQ(ierr);

    } else if ( (solver->PFunction) && ((!solver->JFunction) || (context.flag_jfnk_pre)) ) {

      /* Physical model specifies a preconditioning function, but Jacobian function is not
         specified or user input flag specifies:
         Use: Preconditioned Jacobian-free Newton-Krylov
      */

      if (!mpi->rank) printf("Using the preconditioned Jacobian-free Newtown-Krylov approach.\n");
      /* Preconditioning matrix */
      flag_mat_b = 1;
      ierr = MatCreate  (MPI_COMM_WORLD,&B);                                        CHKERRQ(ierr);
      ierr = MatSetSizes(B,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE);  CHKERRQ(ierr);
      ierr = MatSetBlockSize(B,solver->nvars);                                      CHKERRQ(ierr);
      ierr = MatSetType (B,MATAIJ);                                                 CHKERRQ(ierr);
      ierr = MatSetUp   (B);                                                        CHKERRQ(ierr);
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,B,PetscIJacobianIMEX_JFNK_Pre,&context);           CHKERRQ(ierr);

      SNES snes;
      SNESType snestype;
      ierr = TSGetSNES(ts,&snes);                                                   CHKERRQ(ierr);
      ierr = SNESGetType(snes,&snestype);                                           CHKERRQ(ierr);
      ierr = MatCreateShell(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,
                            PETSC_DETERMINE,&context,&A);                           CHKERRQ(ierr);
      if (!strcmp(snestype,SNESKSPONLY)) {
        context.flag_is_linear = 1;
        ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_Linear);
                                                                                    CHKERRQ(ierr);
      } else {
        context.flag_is_linear = 0;
        ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_JFNK);
                                                                                    CHKERRQ(ierr);
      }
      ierr = MatSetUp(A);                                                           CHKERRQ(ierr);

    } else if ( (solver->PFunction) && (solver->JFunction) ) {

      /* Physical model does specifies Jacobian and preconditioning functions
         Use: the specified Jacobian and preconditioning matrices
      */

      if (!mpi->rank) printf("Using specified Jacobian and preconditioner matrices.\n");
      /* Jacobian matrix */
      ierr = MatCreate  (MPI_COMM_WORLD,&A);                                        CHKERRQ(ierr);
      ierr = MatSetSizes(A,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE);  CHKERRQ(ierr);
      ierr = MatSetBlockSize(A,solver->nvars);                                      CHKERRQ(ierr);
      ierr = MatSetType (A,MATAIJ);                                                 CHKERRQ(ierr);
      ierr = MatSetUp   (A);                                                        CHKERRQ(ierr);
      /* Preconditioning matrix */
      flag_mat_b = 1; 
      ierr = MatCreate  (MPI_COMM_WORLD,&B);                                        CHKERRQ(ierr);
      ierr = MatSetSizes(B,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE);  CHKERRQ(ierr);
      ierr = MatSetBlockSize(B,solver->nvars);                                      CHKERRQ(ierr);
      ierr = MatSetType (B,MATAIJ);                                                 CHKERRQ(ierr);
      ierr = MatSetUp   (B);                                                        CHKERRQ(ierr);
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,B,PetscIJacobianIMEX_Jac_Pre,&context);            CHKERRQ(ierr);

    } else if ( (solver->JFunction) && (context.flag_jfnk_pre) ) {

      /* Physical model does specifies Jacobian function but not preconditioning function, and
         input flag wants to use preconditioned JFNK,
         Use: Preconditioned Jacobian-free Newton-Krylov with the Jacobian function as the 
              preconditioner.
      */

      if (!mpi->rank) printf("Using specified Jacobian as preconditioner to Jacobian-free Newton-Krylov approach.\n");
      /* Preconditioning matrix */
      flag_mat_b = 1; 
      ierr = MatCreate  (MPI_COMM_WORLD,&B);                                        CHKERRQ(ierr);
      ierr = MatSetSizes(B,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE);  CHKERRQ(ierr);
      ierr = MatSetBlockSize(B,solver->nvars);                                      CHKERRQ(ierr);
      ierr = MatSetType (B,MATAIJ);                                                 CHKERRQ(ierr);
      ierr = MatSetUp   (B);                                                        CHKERRQ(ierr);
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,B,PetscIJacobianIMEX_JFNK_JacIsPre,&context);      CHKERRQ(ierr);

      SNES snes;
      SNESType snestype;
      ierr = TSGetSNES(ts,&snes);                                                   CHKERRQ(ierr);
      ierr = SNESGetType(snes,&snestype);                                           CHKERRQ(ierr);
      ierr = MatCreateShell(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,
                            PETSC_DETERMINE,&context,&A);                           CHKERRQ(ierr);
      if (!strcmp(snestype,SNESKSPONLY)) {
        context.flag_is_linear = 1;
        ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_Linear);
                                                                                    CHKERRQ(ierr);
      } else {
        context.flag_is_linear = 0;
        ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_JFNK);
                                                                                    CHKERRQ(ierr);
      }
      ierr = MatSetUp(A);                                                           CHKERRQ(ierr);

    } else {

      /* Physical model does specifies Jacobian function but not preconditioning function
         Use: specified Jacobian matrix, and the same matrix as the preconditioner matrix
      */

      if (!mpi->rank) printf("Using specified Jacobian and using it as a preconditioner too.\n");
      /* Jacobian matrix */
      ierr = MatCreate  (MPI_COMM_WORLD,&A);                                        CHKERRQ(ierr);
      ierr = MatSetSizes(A,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE);  CHKERRQ(ierr);
      ierr = MatSetBlockSize(A,solver->nvars);                                      CHKERRQ(ierr);
      ierr = MatSetType (A,MATAIJ);                                                 CHKERRQ(ierr);
      ierr = MatSetUp   (A);                                                        CHKERRQ(ierr);
      /* Set the IJacobian function for TS */
      ierr = TSSetIJacobian(ts,A,A,PetscIJacobianIMEX_Jac_NoPre,&context);          CHKERRQ(ierr);

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
    ierr = PetscOptionsGetBool(PETSC_NULL,"-parabolic_explicit",&flag,PETSC_NULL);  CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_parabolic = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-parabolic_implicit",&flag,PETSC_NULL);  CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_parabolic = _IMPLICIT_; 

    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-source_explicit",&flag,PETSC_NULL);     CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_source = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-source_implicit",&flag,PETSC_NULL);     CHKERRQ(ierr);
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

  } else {
    fprintf(stderr,"Error in SolvePETSc: TSType %s not supported.\n",time_scheme);
    return(1);
  }

  /* Set pre/post-stage and post-timestep function */
  ierr = TSSetPreStep (ts,PetscPreTimeStep );                             CHKERRQ(ierr);
  ierr = TSSetPreStage(ts,PetscPreStage    );                             CHKERRQ(ierr);
  ierr = TSSetPostStage(ts,PetscPostStage  );                             CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,PetscPostTimeStep);                             CHKERRQ(ierr);
  /* Set solution vector for TS */
  ierr = TSSetSolution(ts,Y);                                             CHKERRQ(ierr);
  /* Set it all up */
  ierr = TSSetUp(ts);                                                     CHKERRQ(ierr);
  /* Set application context */
  ierr = TSSetApplicationContext(ts,&context);                            CHKERRQ(ierr);

  if (!mpi->rank) {
    if (context.flag_is_linear) {
      printf("SolvePETSc(): Note that SNES type is set to %s. ",SNESKSPONLY);
      printf("A linear system will be assumed when (if) computing Jacobian.\n");
    }
  }

  if (!mpi->rank) printf("** Starting PETSc time integration **\n");
  ierr = TSSolve(ts,Y);                                                   CHKERRQ(ierr);
  if (!mpi->rank) printf("** Completed PETSc time integration **\n");

  /* Get the number of time steps */
  ierr = TSGetTimeStepNumber(ts,&solver->n_iter);                         CHKERRQ(ierr);

  /* copy final solution from PETSc's vector */
  ierr = TransferVecFromPETSc(solver->u,Y,&context);                      CHECKERR(ierr);

  /* clean up */
  if (!strcmp(time_scheme,TSARKIMEX)) {
    ierr = MatDestroy(&A);                                                CHKERRQ(ierr);
    if (flag_mat_b) ierr = MatDestroy(&B);                                CHKERRQ(ierr);
  }
  ierr = TSDestroy(&ts);                                                  CHKERRQ(ierr);
  ierr = VecDestroy(&Y);                                                  CHKERRQ(ierr);

  /* write a final solution file, if last iteration did not write one */
  if (context.tic) { IERR OutputSolution(solver,mpi); CHECKERR(ierr); }
  /* calculate error if exact solution has been provided */
  IERR CalculateError(solver,mpi); CHECKERR(ierr);

  PetscFunctionReturn(0);
}

#endif
