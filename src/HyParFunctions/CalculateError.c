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

  size = 1;
  for (i = 0; i < solver->ndims; i++) 
    size *= (solver->dim_local[i]+2*solver->ghosts);


  IERR ExactSolution(solver,mpi,uex,&exact_flag); CHECKERR(ierr);
  if (!exact_flag) {
    /* No exact solution */
    solver->error[0] = solver->error[1] = solver->error[2] = 0.0;
    return(0);  
  }

  /* compute error = difference between exact and numerical solution */
  _ArrayAXPY_(solver->u,-1.0,uex,size*solver->nvars);

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

  free(uex);
  return(0);
}

