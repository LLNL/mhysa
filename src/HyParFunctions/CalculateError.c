#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int CalculateError(void *s,void *m)
{
  HyPar         *solver     = (HyPar*)        s;
  MPIVariables  *mpi        = (MPIVariables*) m;
  int           exact_flag  = 0, d, i, ferr, total_size;
  double        sum         = 0, global_sum = 0;
  _DECLARE_IERR_;

  /* Allocate arrays for reading in the exact solution */
  double *ug,*xg; 
  if (!mpi->rank) {
    int size,offset;

    if (!strcmp(solver->ip_file_type,"ascii")) {

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

        printf("Reading grid and exact solution from ASCII file \"exact.inp\".\n");
        /* read grid, though not necessary (so that file format is consistent) */
        offset = 0;
        for (d = 0; d < solver->ndims; d++) {
          for (i = 0; i < solver->dim_global[d]; i++) ferr = fscanf(in,"%lf",&xg[i+offset]);
          offset += solver->dim_global[d];
        }
        /* read solution */
        for (i = 0; i < solver->nvars; i++) {
          int *index = solver->index;
          int done = 0; _ArraySetValue_(index,solver->ndims,0);
          while (!done) {
            int p; _ArrayIndex1D_(solver->ndims,solver->dim_global,index,0,p);
            ferr = fscanf(in,"%lf",&ug[p*solver->nvars+i]);
            if (ferr != 1) return(1);
            _ArrayIncrementIndex_(solver->ndims,solver->dim_global,index,done);
          }
        }
        fclose(in);

        /* write the exact solution to file in the same format as the computed solution */
        /* to enable comparisons and plotting together                                  */
        if (solver->WriteOutput) {
          char filename[_MAX_STRING_SIZE_] = "op_exact.dat";
          IERR solver->WriteOutput(solver->ndims,solver->nvars,solver->dim_global,xg,ug,
                                   filename,solver->index); CHECKERR(ierr);
        }
      }

    } else if ((!strcmp(solver->ip_file_type,"bin")) || (!strcmp(solver->ip_file_type,"binary"))) {

      FILE *in; in = fopen("exact.inp","rb");
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

        printf("Reading grid and exact solution from binary file \"exact.inp\".\n");

        /* read grid, though not necessary (so that file format is consistent) */
        total_size = 0;
        for (d = 0; d < solver->ndims; d++) total_size += solver->dim_global[d];
        fread(xg, sizeof(double), total_size, in);

        /* read solution */
        total_size = 1;
        for (d = 0; d < solver->ndims; d++) total_size *= solver->dim_global[d]; total_size *= solver->nvars;
        fread(ug, sizeof(double), total_size, in);

        fclose(in);

        /* write the exact solution to file in the same format as the computed solution */
        /* to enable comparisons and plotting together                                  */
        if (solver->WriteOutput) {
          char filename[_MAX_STRING_SIZE_] = "op_exact.bin";
          IERR solver->WriteOutput(solver->ndims,solver->nvars,solver->dim_global,xg,ug,
                                   filename,solver->index); CHECKERR(ierr);
        }
      }

    }

  } else ug = xg = NULL;

  /* Broadcast exact_flag to all processes */
  IERR MPIBroadcast_integer(&exact_flag,1,0,&mpi->world);

  if (!exact_flag)  return(0);  /* No exact solution */

  /* allocate local exact solution array */
  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);
  double *uex = (double*) calloc (size*solver->nvars,sizeof(double));
  
  /* partition global exact solution array to all processes */
  IERR MPIPartitionArraynD(solver->ndims,mpi,(mpi->rank?NULL:ug),uex,
                             solver->dim_global,solver->dim_local,
                             solver->ghosts,solver->nvars); CHECKERR(ierr);

  /* free up global exact solution arrays */
  if (ug)  free(ug);
  if (xg)  free(xg);

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

