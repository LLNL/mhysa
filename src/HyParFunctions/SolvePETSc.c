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

int SolvePETSc(void *s,void *m)
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  PetscErrorCode  ierr    = 0, d;
  TS              ts; /* time integration object  */
  Vec             Y;  /* PETSc solution vector    */
  Mat             A;  /* Jacobian matrix          */
  TSType          time_scheme;

  /* Register custom time-integration methods, if specified */
  ierr = PetscRegisterTIMethods(mpi->rank);                               CHECKERR(ierr);
  if(!mpi->rank) printf("Setting up PETSc time integration... \n");

  /* create and set a PETSc context */
  PETScContext context;
  context.solver = solver;
  context.mpi    = mpi;

  /* create and initialize PETSc solution vector and other parameters */
  /* PETSc solution vector does not have ghost points */
  int total_size = 1;
  for (d=0; d<solver->ndims; d++) total_size *= (solver->dim_local[d]);
  total_size *= solver->nvars;
  ierr = VecCreate(MPI_COMM_WORLD,&Y);                                    CHKERRQ(ierr);
  ierr = VecSetSizes(Y,total_size,PETSC_DECIDE);                          CHKERRQ(ierr);
  ierr = VecSetUp(Y);                                                     CHKERRQ(ierr);

  /* copy initial solution to PETSc's vector */
  ierr = TransferToPETSc(solver->u,Y,&context);                           CHECKERR(ierr);

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
    ierr = MatCreateShell(MPI_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,
                          PETSC_DETERMINE,&context,&A);                   CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX);
                                                                          CHKERRQ(ierr);
    ierr = MatSetUp(A);                                                   CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts,A,A,PetscIJacobianIMEX,&context);            CHKERRQ(ierr);

    /* set pre-conditioner to none for MatShell */
    SNES snes;
    KSP  ksp;
    PC   pc;
    TSGetSNES(ts,&snes);
    SNESGetKSP(snes,&ksp);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCNONE);

    /* read the implicit/explicit flags for each of the terms for IMEX schemes */
    /* default -> hyperbolic - explicit, parabolic and source - implicit       */
    PetscBool flag = PETSC_FALSE;

    context.flag_hyperbolic = _EXPLICIT_; 
    context.flag_parabolic  = _IMPLICIT_; 
    context.flag_source     = _IMPLICIT_; 

    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_hyperbolic = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_hyperbolic = _IMPLICIT_; 

    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic1_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_hyperbolic1 = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic1_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_hyperbolic1 = _IMPLICIT_; 

    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic2_explicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_hyperbolic2 = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    ierr = PetscOptionsGetBool(PETSC_NULL,"-hyperbolic2_implicit",&flag,PETSC_NULL); CHKERRQ(ierr);
    if (flag == PETSC_TRUE) context.flag_hyperbolic2 = _IMPLICIT_; 

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

  } else {
    fprintf(stderr,"Error in SolvePETSc: TSType %s not supported.\n",time_scheme);
    return(1);
  }

  /* Set pre/post-stage and post-timestep function */
  ierr = TSSetPreStage(ts,PetscPreStage    );                             CHKERRQ(ierr);
  ierr = TSSetPostStage(ts,PetscPostStage  );                             CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,PetscPostTimeStep);                             CHKERRQ(ierr);
  /* Set solution vector for TS */
  ierr = TSSetSolution(ts,Y);                                             CHKERRQ(ierr);
  /* Set it all up */
  ierr = TSSetUp(ts);                                                     CHKERRQ(ierr);
  /* Set application context */
  ierr = TSSetApplicationContext(ts,&context);                            CHKERRQ(ierr);

  if (!mpi->rank) printf("** Starting PETSc time integration **\n");
  ierr = TSSolve(ts,Y);                                                   CHKERRQ(ierr);
  if (!mpi->rank) printf("** Completed PETSc time integration **\n");

  /* copy final solution from PETSc's vector */
  ierr = TransferFromPETSc(solver->u,Y,&context);

  /* clean up */
  if (!strcmp(time_scheme,TSARKIMEX)) {
    ierr = MatDestroy(&A);                                                CHKERRQ(ierr);
  }
  ierr = TSDestroy(&ts);                                                  CHKERRQ(ierr);
  ierr = VecDestroy(&Y);                                                  CHKERRQ(ierr);
  return(0);
}

#endif
