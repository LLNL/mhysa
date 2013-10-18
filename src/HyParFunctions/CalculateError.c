#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int CalculateError(void *s,void *m)
{
  HyPar         *solver     = (HyPar*)        s;
  MPIVariables  *mpi        = (MPIVariables*) m;
  int           exact_flag  = 0, d, i, ierr = 0;
  double        sum         = 0, global_sum = 0;

  /* Allocate arrays for reading in the exact solution */
  double *ug,*xg; 
  if (!mpi->rank) {
    int size,offset;

    /* Reading grid and exact solution */
    FILE *in; in = fopen("exact.inp","r");
    if (!in) {
      solver->error[0] = solver->error[1] = solver->error[2] = 0;
      exact_flag = 0;
      ug = xg = NULL;
    } else {
      exact_flag = 1;
      /* allocate global solution array */
      size = solver->npoints_global*solver->nvars;
      ug = (double*) calloc(size,sizeof(double));
      size = 0;
      for (d=0; d<solver->ndims; d++) size += solver->dim_global[d];
      xg      = (double*) calloc(size,sizeof(double));

      printf("Reading grid and exact solution from file \"exact.inp\".\n");
      /* read grid, though not necessary (so that file format is consistent) */
      offset = 0;
      for (d = 0; d < solver->ndims; d++) {
        for (i = 0; i < solver->dim_global[d]; i++) ierr = fscanf(in,"%lf",&xg[i+offset]);
        offset += solver->dim_global[d];
      }
      /* read solution */
      for (i = 0; i < solver->nvars; i++) {
        int *index = solver->index;
        int done = 0; ierr = ArraySetValue_int(index,solver->ndims,0); CHECKERR(ierr);
        while (!done) {
          int p = ArrayIndex1D(solver->ndims,solver->dim_global,index,NULL,0);
          ierr = fscanf(in,"%lf",&ug[p*solver->nvars+i]);
          if (ierr != 1) return(1);
          done = ArrayIncrementIndex(solver->ndims,solver->dim_global,index);
        }
      }
      fclose(in);

      /* write the exact solution to file in the same format as the computed solution */
      /* to enable comparisons and plotting together                                  */
      char filename[_MAX_STRING_SIZE_] = "op_exact.dat";
      ierr = solver->WriteOutput(solver->ndims,solver->nvars,solver->dim_global,xg,ug,
                                 filename,solver->index); CHECKERR(ierr);
    }

  } else ug = xg = NULL;

  /* Broadcast exact_flag to all processes */
  ierr = MPIBroadcast_integer(&exact_flag,1,0,&mpi->world);

  if (!exact_flag)  return(0);  /* No exact solution */

  /* allocate local exact solution array */
  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);
  double *uex = (double*) calloc (size*solver->nvars,sizeof(double));
  
  /* partition global exact solution array to all processes */
  ierr = MPIPartitionArraynD(solver->ndims,mpi,(mpi->rank?NULL:ug),uex,
                             solver->dim_global,solver->dim_local,
                             solver->ghosts,solver->nvars); CHECKERR(ierr);

  /* free up global exact solution arrays */
  if (ug)  free(ug);
  if (xg)  free(xg);

  /* compute error = difference between exact and numerical solution */
  ierr = ArrayAXPY(solver->u,-1.0,uex,size*solver->nvars); CHECKERR(ierr);

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

