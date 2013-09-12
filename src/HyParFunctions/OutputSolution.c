#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/* Function declarations */
static int IncrementFilename (char*);

int OutputSolution(void *s, void *m)
{
  HyPar         *solver = (HyPar*)       s;
  MPIVariables  *mpi    = (MPIVariables*)m;
  int           ierr    = 0,d;

  /* if WriteOutput() is NULL, then return */
  if (!solver->WriteOutput) return(0);

  /* root process: allocate global output arrays */
  double *ug, *xg;
  if (!mpi->rank) {
    int size_global;
    size_global = 1;
    for (d=0; d<solver->ndims; d++) size_global *= solver->dim_global[d];
    ug = (double*) calloc (size_global,sizeof(double));
    ierr = ArraySetValue_double(ug,size_global,0.0); CHECKERR(ierr);
    size_global = 0;
    for (d=0; d<solver->ndims; d++) size_global += solver->dim_global[d];
    xg = (double*) calloc (size_global,sizeof(double));
    ierr = ArraySetValue_double(xg,size_global,0.0); CHECKERR(ierr);
  } else {
    /* null pointers on non-root processes */
    ug = xg = NULL;
  }

  /* Assemble the local output arrays into the global output arrays */
  ierr = MPIGatherArraynD(solver->ndims,mpi,ug,solver->u,solver->dim_global,
                          solver->dim_local,solver->ghosts,solver->nvars);  CHECKERR(ierr);
  int offset_global, offset_local;
  offset_global = offset_local = 0;
  for (d=0; d<solver->ndims; d++) {
    ierr = MPIGatherArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),&solver->x[offset_local],
                            mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    offset_global += solver->dim_global[d];
    offset_local  += solver->dim_local [d];
  }

  if (!mpi->rank) {
    /* write output file to disk */
    char filename[_MAX_STRING_SIZE_] = "";
    strcat(filename,solver->op_filename);
    ierr = solver->WriteOutput(solver->ndims,solver->nvars,solver->dim_global,xg,ug,
                               filename,solver->index); CHECKERR(ierr);
    if (!strcmp(solver->op_overwrite,"no"))  ierr = IncrementFilename(solver->op_filename);

    /* Clean up output arrays */
    free(xg);
    free(ug);
  }

  return(0);
}

int IncrementFilename(char *f)
{
  int ierr = 0;
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
            ierr = 1;
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
  return(ierr);
}
