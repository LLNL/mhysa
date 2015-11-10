/*! @file CalculateError.c
    @author Debojyoti Ghosh
    @brief Computes the error in the solution.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <timeintegration.h>
#include <mpivars.h>
#include <hypar.h>

int ExactSolution(void*,void*,double*,int*);

/*! Calculates the error in the solution if the exact solution is 
    available. If the exact solution is not available, the errors
    are reported as zero.
    The exact solution should be provided in the file "exact.inp"
    in the same format as the initial solution.
*/
int CalculateError(
                    void *s, /*!< Solver object of type #HyPar */
                    void *m  /*!< MPI object of type #MPIVariables */
                  )
{
  HyPar         *solver     = (HyPar*)        s;
  MPIVariables  *mpi        = (MPIVariables*) m;
  int           exact_flag  = 0, i, size;
  double        sum         = 0, global_sum = 0;
  double        *uex        = NULL;
  _DECLARE_IERR_;

  size = solver->nvars;
  for (i = 0; i < solver->ndims; i++) 
    size *= (solver->dim_local[i]+2*solver->ghosts);
  uex = (double*) calloc (size, sizeof(double));

  static const double tolerance = 1e-15;
  IERR ExactSolution(solver,mpi,uex,&exact_flag); CHECKERR(ierr);

  if (!exact_flag) {

    /* No exact solution */
    IERR TimeError(solver,mpi,NULL); CHECKERR(ierr);
    solver->error[0] = solver->error[1] = solver->error[2] = 0.0;

  } else {

    IERR TimeError(solver,mpi,uex); CHECKERR(ierr);

    /* calculate solution norms (for relative error) */
    double solution_norm[3] = {0.0,0.0,0.0};
    /* L1 */
    sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,uex);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    solution_norm[0] = global_sum/((double)solver->npoints_global);
    /* L2 */
    sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,uex);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    solution_norm[1] = sqrt(global_sum/((double)solver->npoints_global));
    /* Linf */
    sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,uex);
    global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
    solution_norm[2] = global_sum;

    /* compute error = difference between exact and numerical solution */
    _ArrayAXPY_(solver->u,-1.0,uex,size);

    /* calculate L1 norm of error */
    sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,uex);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    solver->error[0] = global_sum/((double)solver->npoints_global);

    /* calculate L2 norm of error */
    sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,uex);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    solver->error[1] = sqrt(global_sum/((double)solver->npoints_global));

    /* calculate Linf norm of error */
    sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,uex);
    global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
    solver->error[2] = global_sum;

    /* 
      decide whether to normalize and report relative errors, 
      or report absolute errors.
    */
    if (    (solution_norm[0] > tolerance) 
        &&  (solution_norm[1] > tolerance) 
        &&  (solution_norm[2] > tolerance) ) {
      solver->error[0] /= solution_norm[0];
      solver->error[1] /= solution_norm[1];
      solver->error[2] /= solution_norm[2];
    }
  }

  free(uex);
  return(0);
}

