#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ExactSolution(void*,void*,double*,int*);

int CalculateError(void *s,void *m)
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
    solver->error[0] = solver->error[1] = solver->error[2] = 0.0;

  } else {

    /* calculate solution norms (for relative error) */
    double solution_norm[3] = {0.0,0.0,0.0};
    /* L1 */
    sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,solver->u);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    solution_norm[0] = global_sum/((double)solver->npoints_global);
    /* L2 */
    sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,solver->u);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    solution_norm[1] = sqrt(global_sum/((double)solver->npoints_global));
    /* Linf */
    sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                           solver->ghosts,solver->index,solver->u);
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

