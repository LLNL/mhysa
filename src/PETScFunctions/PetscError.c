#ifdef with_petsc

/*! @file PetscError.c
    @author Debojyoti Ghosh
    @brief Compute time integration error estimates, if available
*/

#include <math.h>
#include <string.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscTimeError"

/*! Compute the norms of the error estimate, if the PETSc time integrator 
    has it (for example TSGLEE)
*/
int PetscTimeError(
                    TS  ts /*!< Time integrator object of PETSc type TS */
                  )
{
  PETScContext    *context  = NULL;
  HyPar           *solver   = NULL;
  MPIVariables    *mpi      = NULL;
  Vec             Y,Z;
  TSType          time_scheme;
  double          dt;
  PetscErrorCode  ierr;

  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscError: Null context!\n");
    return(1);
  }
  solver  = context->solver;
  mpi     = context->mpi;

  int    size    = solver->npoints_local_wghosts * solver->nvars;
  double sum     = 0.0, global_sum = 0.0, *Uerr = solver->uref,
         error[3] = {0,0,0};

  ierr = TSGetTimeStep(ts,&dt  ); CHKERRQ(ierr);
  ierr = TSGetType(ts,&time_scheme); CHKERRQ(ierr);
  ierr = TSGetSolution(ts,&Y); CHKERRQ(ierr);
  ierr = TransferVecFromPETSc(solver->u,Y,context); CHECKERR(ierr);
  if (!strcmp(time_scheme,TSGLEE)) {

    ierr = VecDuplicate(Y,&Z); CHKERRQ(ierr);
    ierr = TSGetTimeError(ts,0,&Z);CHKERRQ(ierr);
    ierr = TransferVecFromPETSc(Uerr,Z,context); CHECKERR(ierr);
    ierr = VecDestroy(&Z); CHKERRQ(ierr);

    /* calculate solution norm for relative errors */
    double sol_norm[3] = {0.0,0.0,0.0};
    /* L1 */
    sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,solver->u);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    sol_norm[0] = global_sum/((double)solver->npoints_global);
    /* L2 */
    sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,solver->u);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    sol_norm[1] = sqrt(global_sum/((double)solver->npoints_global));
    /* Linf */
    sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,solver->u);
    global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
    sol_norm[2] = global_sum;

    /* calculate L1 norm of error */
    sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,Uerr);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    error[0] = global_sum/((double)solver->npoints_global);
    /* calculate L2 norm of error */
    sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,Uerr);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    error[1] = sqrt(global_sum/((double)solver->npoints_global));
    /* calculate Linf norm of error */
    sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,Uerr);
    global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
    error[2] = global_sum;

    if (   (sol_norm[0] > _MACHINE_ZERO_)
        && (sol_norm[1] > _MACHINE_ZERO_)
        && (sol_norm[2] > _MACHINE_ZERO_) ) {
      error[0] /= sol_norm[0];
      error[1] /= sol_norm[1];
      error[2] /= sol_norm[2];
    }

    /* write to file */
    if (!mpi->rank) {
      FILE *out;
      out = fopen("glm_err.dat","w");
      fprintf(out,"%1.16E  %1.16E  %1.16E  %1.16E  ",dt,error[0],error[1],error[2]);
      fclose(out);
      printf("Estimated time integration errors (GLM-GEE time-integration):-\n");
      printf("  L1         Error           : %1.16E\n",error[0]);
      printf("  L2         Error           : %1.16E\n",error[1]);
      printf("  Linfinity  Error           : %1.16E\n",error[2]);
    }
  }

  return(0);
}

#endif
